#include <iostream>
#include <fstream>
#include <cstring>

#include "bwtgap.h"
#include "bwtaln.h"
#include "bwamem.h"

#include "CTPL/ctpl_stl.h"

#include "Fastq.hpp"
#include "QualityScores.hpp"

//------------------ Data structures ----------------

typedef std::vector<SGLib::Fastq::FastqRecord> yalff_chunk_t;
typedef yalff_chunk_t::iterator chunk_itr;
typedef std::pair<std::size_t, std::size_t> chunk_view;

struct yalff_opt_t{
    std::size_t k; //k-mer length
    std::size_t mismatches; //number of mismatches allowed
    std::size_t n_threads; //number of threads
    std::size_t chunk_size; //number of fastq records read at each round
    char bad_threshold;
    char good_threshold;
    char replacement; //smoothing character
    char bad_replacement;
    std::size_t skip; //number of bases to skip
    bwaidx_t* index;
};

class semaphore{

    public:
        semaphore(std::size_t resources);
        void wait(std::size_t to_take);
        void signal(std::size_t to_give);
        ~semaphore();

    private:
        std::mutex mtx;
        std::condition_variable cond;
        std::size_t counter;
        const std::size_t global_resources;
};

typedef uint8_t k_t;
const k_t k_max = std::numeric_limits<k_t>::max(); //not included

//------------------ Prototypes ---------------------

yalff_opt_t yalff_init(int argc, char* argv[]);

void yalff_destroy(const yalff_opt_t& yopt);

inline gap_opt_t getAlnOpts(const yalff_opt_t& yopt);

yalff_chunk_t read_chunk(std::istream& in, const yalff_opt_t& yopt);

int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);

std::pair<bwt_aln1_t*, int> getKmerMatches(bwt_t *bwt, int qlen, uint8_t* query, const gap_opt_t& opt);

void smooth_chunk(int thread_id, const yalff_opt_t& yopt, semaphore& s, yalff_chunk_t& slice, std::size_t start, std::size_t stop);

void countBadQuals(const std::string& qs, std::size_t k, char thr, std::vector<k_t>& toRet);

//---------------------------------------------------

int main(int argc, char* argv[])
{
    //std::cerr << "Begin program" << std::endl;

    const yalff_opt_t yopt = yalff_init(argc, argv);
    //std::cerr << "Options read" << std::endl;
    if(yopt.k == 0) return 1;//error

    ctpl::thread_pool pool(static_cast<int>(yopt.n_threads)); //thread pool
    semaphore sema(yopt.n_threads);
    uint8_t chunk_counter = 0;
    std::vector<yalff_chunk_t> chunks(2);//rotating buffer

    bool is_first = true;

    while(std::cin.good())
    {
        //std::cerr << "Reading chunk" << std::endl;

        //read chunk
        chunks[chunk_counter] = read_chunk(std::cin, yopt);

        //std::cerr << "Chunk read" << std::endl;

        //subdivide it into subgroups
        std::vector<chunk_view> intervals(yopt.n_threads);
        std::size_t step_index = chunks[chunk_counter].size() % yopt.n_threads;
        std::size_t sch_size = chunks[chunk_counter].size() / yopt.n_threads + 1;
        for(std::size_t i = 0; i < step_index; ++i)//if step_index == 0 then chunk_size is divisible by n_threads -> skip this loop
        {
            intervals.push_back(chunk_view(i * sch_size, (i + 1) * sch_size));
        }
        --sch_size; //the second part of the views have size smaller than one
        for(std::size_t i = 0; i < yopt.n_threads - step_index; ++i)
        {
            intervals.push_back(chunk_view(step_index * (sch_size + 1) + i * sch_size, step_index * (sch_size + 1) + (i + 1) * sch_size));
        }

        //std::cerr << "step index = " << step_index << std::endl;
        //std::cerr << "sch size = " << sch_size << std::endl;
        //std::cerr << "Chunk subdivided" << std::endl;

        //wait until the previous chunk has been processed
        sema.wait(yopt.n_threads);
        //std::cerr << "Previous chunk processed" << std::endl;
        sema.signal(yopt.n_threads);

        //use the thread pool to smooth it in parallel while writing the previous chunk
        for(auto interval : intervals)
        {
            pool.push(smooth_chunk, std::ref(yopt), std::ref(sema), std::ref(chunks[chunk_counter]), interval.first, interval.second);
        }

        //std::cerr << "Chunk processed" << std::endl;

        if(chunk_counter == 0) chunk_counter = 1; //select the other chunk (already processed)
        else chunk_counter = 0;

        //std::cerr << "Writing chunk" << std::endl;

        for(auto& record : chunks[chunk_counter])
        {
            if(is_first)
            {
                is_first = false;
                std::cout << record;
            }
            else std::cout << "\n" << record;
        }

        //std::cerr << "Chunk written" << std::endl;

    }

    sema.wait(yopt.n_threads);
    //std::cerr << "Last chunk processed" << std::endl;
    sema.signal(yopt.n_threads);
    //std::cerr << "Stopping pool" << std::endl;
    pool.stop(true);//wait for the current threads to stop (the queue is cleared regardless to the waiting jobs)

    //save the last chunk
    if(chunk_counter == 0) chunk_counter = 1;
    else chunk_counter = 0;
    for(auto& record : chunks[chunk_counter])
    {
        if(is_first)
        {
            is_first = false;
            std::cout << record;
            //std::cerr << record;
        }
        else
        {
            std::cout << "\n" << record;
            //std::cerr << "\n" << record;
        }
    }

    //std::cerr << std::endl;

    yalff_destroy(yopt);
    return 0;
}

//------------------ Implementations ---------------

semaphore::semaphore(std::size_t resources) : global_resources(resources)
{
    counter = resources;
}

void semaphore::wait(std::size_t to_take)
{
    std::unique_lock<decltype(mtx)> lock(mtx);
    while(counter < to_take) cond.wait(lock);
    counter -= to_take;
}

void semaphore::signal(std::size_t to_give)
{
    std::unique_lock<decltype(mtx)> lock(mtx);
    counter += to_give;
    if(counter > global_resources) counter = global_resources;
    cond.notify_all();
}

semaphore::~semaphore()
{
    //dtor
}

//reminder: Sanger scale = !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
yalff_opt_t yalff_init(int argc, char* argv[])
{
    yalff_opt_t opts;
    //default options
    opts.chunk_size = 10000;
    opts.k = 32;
    opts.mismatches = 1;
    opts.n_threads = std::thread::hardware_concurrency() != 1 ? std::thread::hardware_concurrency() - 1 : 1;
    opts.replacement = 'I';
    opts.bad_replacement = 'j';
    opts.bad_threshold = '$';
    opts.good_threshold = 'I';
    opts.skip = 0;

    //read the parameters
    for(int i = 1; i < argc; ++i)
    {
        if(strcmp(argv[i], "-k") == 0)
        {
            //k-mer length
            ++i;
            opts.k = strtoull(argv[i], nullptr, 10);
        }
        else if(strcmp(argv[i], "-m") == 0)
        {
            //number of mismatches allowed in a k-mer
            ++i;
            opts.mismatches = strtoull(argv[i], nullptr, 10);
        }
        else if(strcmp(argv[i], "-c") == 0)
        {
            //chunk size
            ++i;
            opts.chunk_size = strtoull(argv[i], nullptr, 10);
        }
        else if(strcmp(argv[i], "-q") == 0)
        {
            //smoothing value
            ++i;
            opts.replacement = *(argv[i]);
        }
        else if(strcmp(argv[i], "-e") == 0)
        {
            //smoothing value (bad)
            ++i;
            opts.bad_replacement = *(argv[i]);
        }
        else if(strcmp(argv[i], "-b") == 0)
        {
            //bad threshold
            ++i;
            opts.bad_threshold = *(argv[i]);
        }
        else if(strcmp(argv[i], "-g") == 0)
        {
            //good threshold
            ++i;
            opts.good_threshold = *(argv[i]);
        }
        else if(strcmp(argv[i], "-s") == 0)
        {
            //number of bases to skip
            ++i;
            opts.skip = strtoull(argv[i], nullptr, 10);
        }
        else if(strcmp(argv[i], "-d") == 0)
        {
            //indexed dictionary string
            ++i;
            opts.index = bwa_idx_load(argv[i], BWA_IDX_ALL);
        }
        else if(strcmp(argv[i], "-shm") == 0)
        {
            //indexed dictionary string from shared memory
            ++i;
            opts.index = bwa_idx_load_from_shm(argv[i]);
        }
        else if(strcmp(argv[i], "-t") == 0)
        {
            //number of threads available
            ++i;
            opts.n_threads = strtoull(argv[i], nullptr, 10);
            if(opts.n_threads == 0) opts.n_threads = 1;
        }
        else // -h option
        {
            //usage
            std::cerr << "OPTIONS:\n\n";
            std::cerr << "-k NUM\t k-mer length. [32] (max = " << static_cast<uint32_t>(k_max) << " excluded)\n\n";
            std::cerr << "-m NUM\t Number of mismatches allowed per k-mer. [1]\n\n";
            std::cerr << "-c NUM\t Chunk size. The number of reads read at once on each iteration. [10000]\n\n";
            std::cerr << "-b CHAR\t Sanger threshold for a quality score to be considered. [$]\n\n";
            std::cerr << "-g CHAR\t Sanger threshold for a quality score to be considered correct independently from the dictionary. [I]\n\n";
            std::cerr << "-s NUM\t Number of bases to skip after each k-mer. A value of 0 checks all the k-mers. [0]\n\n";
            std::cerr << "-q CHAR\t Sanger value used as replacement during smoothing. [I]\n\n";
            std::cerr << "-e CHAR\t Sanger value used as an eventual replacement when a k-mer aligns badly. [j]\n\n";
            std::cerr << "-d STRING Path to the indexed fasta (same format used for BWA).\n\n";
            std::cerr << "-shm STRING\t Path to the indexed fasta (same format used for BWA).\n";
            std::cerr << "\t\t The FM-Index MUST be already present shared memory.\n";
            std::cerr << "\t\t It can be loaded and managed using BWA shm command.\n\n";
            std::cerr << "-t NUM\t Number of threads available. [Hardware concurrency - 1]\n\n";
            std::cerr << "-h \t See help.\n\n";
            opts.k = 0;
            i = argc;
        }
    }
    return opts;
}

void yalff_destroy(const yalff_opt_t& yopt)
{
    //free the index if it was loaded by this process
    if(yopt.index->is_shm != 1) bwa_idx_destroy(yopt.index);
}

inline gap_opt_t getAlnOpts(const yalff_opt_t& yopt)
{
    gap_opt_t o;
    /* IMPORTANT: s_mm*10 should be about the average base error rate. Violating this requirement will break pairing! */
    o.s_mm = 3;
    o.s_gapo = 11;
    o.s_gape = 4;
    o.max_diff = 1;
    o.max_gapo = 0;
    o.max_gape = 0;
    o.indel_end_skip = 5;
    o.max_del_occ = 10;
    o.max_entries = std::numeric_limits<int>::max();
    o.mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
    o.seed_len = static_cast<int>(yopt.k);
    o.max_seed_diff = static_cast<int>(yopt.mismatches);
    o.fnr = -1;
    o.n_threads = 1;
    o.max_top2 = 30;
    o.trim_qual = 0;
    if (o.max_diff < o.max_gapo) o.max_gapo = o.max_diff;
    return o;
}

yalff_chunk_t read_chunk(std::istream& in, const yalff_opt_t& yopt)
{
    yalff_chunk_t toRet;
    toRet.reserve(yopt.chunk_size);
    try{
        for(std::size_t i = 0; i < yopt.chunk_size; ++i)
        {
            toRet.emplace_back(SGLib::Fastq::readNextFastqRecord(in));
        }
    } catch(SGLib::fastx_end& file_ended){
        static_cast<void>(file_ended);
    }
    return toRet;
}

// width must be filled as zero
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = bwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = bwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

std::pair<bwt_aln1_t*, int> getKmerMatches(bwt_t *bwt, int qlen, uint8_t* query, const gap_opt_t& opt)
{
    std::pair<bwt_aln1_t*, int> alignments;
    gap_stack_t* stk = gap_init_stack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt);
    bwt_width_t width[qlen + 1];
    bwt_cal_width(bwt, qlen, query, width);
    alignments.first = bwt_match_gap(bwt, qlen, query, width, 0, &opt, &alignments.second, stk);
    gap_destroy_stack(stk);
    return alignments;
}

void smooth_chunk(int thread_id, const yalff_opt_t& yopt, semaphore& s, yalff_chunk_t& slice, std::size_t start, std::size_t stop)
{
    //smooth thread
    s.wait(1);
    for(; start != stop; ++start) //smooth the read
    {
        SGLib::Fastq::FastqRecord& record = slice[start];

        if(record.getSequenceRef().size() >= yopt.k)
        {
            std::vector<uint8_t> read(record.getSequenceRef().size());
            for(std::size_t i = 0; i < record.getSequenceRef().size(); ++i) read[i] = SGLib::base_to_int[static_cast<std::size_t>(record.getSequenceRef()[i])];

            std::vector<bool> to_boost(record.getQualityRef().size(), true); //better cache utilization

            std::vector<k_t> qmer_bad_count(read.size() - yopt.k + 1, 0); //number of bad qualities in each k-mer
            countBadQuals(record.getQualityRef(), yopt.k, yopt.bad_threshold, qmer_bad_count);
            std::vector<uint8_t> qmer_not_so_bad_count(qmer_bad_count.size(), 0);
            countBadQuals(record.getQualityRef(), yopt.k, yopt.good_threshold, qmer_not_so_bad_count);

            for(std::size_t i = 0; i < read.size() - yopt.k + 1; ++i)
            {
                if(qmer_not_so_bad_count[i] != 0 && qmer_bad_count[i] == 0)
                {
                    std::pair<bwt_aln1_t*, int> aln = getKmerMatches(yopt.index->bwt, static_cast<int>(yopt.k), &read[i], getAlnOpts(yopt));
                    bwa_seq_t alignment;
                    bwa_aln2seq(aln.second, aln.first, &alignment); //get sa value
                    if(alignment.type == BWA_TYPE_UNIQUE || alignment.type == BWA_TYPE_REPEAT) //Try for a match with at most yopt.mismatches mismatches
                    {
                        bwtint_t position = bwt_sa(yopt.index->bwt, alignment.sa);
                        int64_t len;
                        uint8_t* ref_aln = bns_get_seq(yopt.index->bns->l_pac, yopt.index->pac, static_cast<int64_t>(position), static_cast<int64_t>(position + yopt.k), &len);
                        if(yopt.k == static_cast<std::size_t>(len)) // possible if out of range or crossed boundary -> error, skip this iteration
                        {
                            for(std::size_t j = 0; j < yopt.k; ++j)//check where the mismatches are and set the to_boost vector accordingly
                            {
                                if(to_boost[i + j] && read[i + j] != ref_aln[j]) to_boost[i + j] = false;
                            }
                        }
                        free(ref_aln);
                    } else qmer_bad_count[i] = k_max;
                    free(aln.first);
                    i += yopt.skip;
                }
            }

            //Boost the quality values based on the qmer_bad_count and to_boost arrays.
            bool good_block_started = false;
            bool bad_block_started = false;
            for(std::size_t i = 0; i < qmer_bad_count.size(); ++i)
            {
                if(qmer_bad_count[i] == 0) //start block of good k-mers
                {
                    if(!good_block_started)
                    {
                        for(std::size_t j = 0; j < yopt.k; ++j) if(to_boost[i + j]) record.setQualityScore(i + j, yopt.replacement);
                        good_block_started = true;
                        bad_block_started = false;
                    }
                    else
                    {
                        if(to_boost[i + yopt.k - 1]) record.setQualityScore(i + yopt.k - 1, yopt.replacement);
                    }
                }
                else if(qmer_bad_count[i] == k_max) //start block of bad k-mers
                {
                    if(!bad_block_started)
                    {
                        for(std::size_t j = 0; j < yopt.k; ++j) if(yopt.bad_replacement < record.getQualityScore(i + j) && to_boost[i + j]) record.setQualityScore(i + j, yopt.bad_replacement);
                        bad_block_started = true;
                        good_block_started = false;
                    }
                    else
                    {
                        if(yopt.bad_replacement < record.getQualityScore(i + yopt.k - 1) && to_boost[i + yopt.k - 1]) record.setQualityScore(i + yopt.k - 1, yopt.bad_replacement);
                    }
                }
                else //close block of good k-mers
                {
                    good_block_started = false;
                    bad_block_started = false;
                }
            }
        }
    }
    s.signal(1);
    static_cast<void>(thread_id);
}

void countBadQuals(const std::string& qs, std::size_t k, char thr, std::vector<k_t>& toRet)
{
    for(std::size_t i = 0; i < k; ++i)
    {
        if(qs[i] < thr)
        {
            ++toRet[0]; //increment count
        }
    }
    for(uint8_t i = 1; i < static_cast<uint8_t>(qs.length() - k + 1); ++i)
    {
        if(qs[i-1] < thr) toRet[i] = toRet[i-1] - 1;
        else toRet[i] = toRet[i-1];
        if(qs[i+k-1] < thr)
        {
            ++toRet[i];
        }
    }
}
