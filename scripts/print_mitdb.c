/**
 * This simple set of utilities prints to stdout the k-mers inside the dictionary of
 * 1) Quartz
 *
 * 2) LAVA
 *    |--> The dictionary of the reference AND
 *    |--> The dictionary containing the k-mers for the SNPs
 *
 * The parameters are -q for Quartz files, -r for LAVA reference, -s for LAVA SNPs k-mer reference.
 * It is possible to specify multiple references in a single command.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define POW_2_32 (1UL << 32)

typedef u_int64_t pac_t;
typedef u_int32_t count_t;
typedef u_int8_t flag_t;

pac_t read_pac(FILE *in)
{
	pac_t x;
	assert(fread(&x, sizeof(pac_t), 1, in) == 1);
	return x;
}

count_t read_count(FILE *in)
{
	count_t x;
	assert(fread(&x, sizeof(count_t), 1, in) == 1);
	return x;
}

flag_t read_flag(FILE *in)
{
	flag_t x;
	assert(fread(&x, sizeof(flag_t), 1, in) == 1);
	return x;
}

void decode_kmer(pac_t e, char* read_array)
{
    int i;
	for (i=0; i<32; ++i)
	{
		switch (e & 0x0000000000000003)
        {
			case 0:
				read_array[i]='A';
				break;
			case 1:
				read_array[i]='C';
				break;
			case 2:
				read_array[i]='G';
				break;
			case 3:
				read_array[i]='T';
				break;
		}
		e >>= 2;
	}
}

int print_quartz_dict(char* path_to_quartz_dict, pac_t* counter)
{
    pac_t packed;
    pac_t dict_siz;
    pac_t i;
    char kmer[33];
    FILE* fp;
    kmer[32] = '\0';
    
    fp = fopen(path_to_quartz_dict, "rb");
    if(!fp) 
    {
        fprintf(stderr, "Unable to open the specified file %s\n", path_to_quartz_dict);
        return 1;
    }   
    
    if(fread((unsigned char*) &dict_siz, 8, 1, fp) != 1)
    {
        fprintf(stderr, "Unable to read the dictionary file %s\n", path_to_quartz_dict);
        return 2;
    }
    
    for(i = 0; i < dict_siz; ++i) {
        fread(&packed, sizeof(pac_t), 1, fp);
        decode_kmer(packed, kmer);
        printf(">k\n%s\n", kmer);
    }
    
    fclose(fp);
    return 0;
}

int print_reference(char* path_to_lava_ref, pac_t* counter)
{
    pac_t dict_siz, aux_table_size;
    pac_t packed;
    count_t kmer_count;
    flag_t amb_flag;
    pac_t i;
    char kmer[33];
    FILE* fp;
    kmer[32] = '\0';
    
    fp = fopen(path_to_lava_ref, "rb");
    if(!fp) 
    {
        fprintf(stderr, "Unable to open the specified ref file %s\n", path_to_lava_ref);
        return 1;
    }   
    
    if(fread((unsigned char*) &dict_siz, sizeof(pac_t), 1, fp) != 1) //read the two uint64_t at the beginning (size and placeholder)
    {
        fprintf(stderr, "Unable to read the total number of kmers for %s\n", path_to_lava_ref);
        return 2;
    }
    if(fread((unsigned char*) &aux_table_size, sizeof(pac_t), 1, fp) != 1)
    {
        fprintf(stderr, "Unable to read the size of the auxiliary table for %s\n", path_to_lava_ref);
        return 3;
    }
    if (dict_siz > POW_2_32) {
		fprintf(stderr, "Reference dictionary is too large (limit: %lu 32-mers)\n", POW_2_32);
		return 4;
    }
    
    for(i = 0; i < dict_siz; ++i) {
        packed = read_pac(fp);
        kmer_count = read_count(fp); //kmer counter
        amb_flag = read_flag(fp); //AMBIGUOUS FLAG
        decode_kmer(packed, kmer);
        printf(">k\n%s\n", kmer); //print the kmer only
    }
    
    fclose(fp);
    return 0;
}

int print_snp_dict(char* path_to_snp_dict, pac_t* counter)
{
    pac_t dict_siz, aux_table_size;
    pac_t packed;
    count_t kmer_count;
    flag_t buffer;
    
    pac_t i;
    char kmer[33];
    FILE* fp;
    kmer[32] = '\0';
    
    fp = fopen(path_to_snp_dict, "rb");
    if(!fp) 
    {
        fprintf(stderr, "Unable to open the specified SNP file %s\n", path_to_snp_dict);
        return 1;
    }   
    
    if(fread((unsigned char*) &dict_siz, sizeof(pac_t), 1, fp) != 1) //read the two uint64_t at the beginning (size and placeholder)
    {
        fprintf(stderr, "Unable to read the total number of kmers for %s\n", path_to_snp_dict);
        return 2;
    }
    if(fread((unsigned char*) &aux_table_size, sizeof(pac_t), 1, fp) != 1)
    {
        fprintf(stderr, "Unable to read the size of the auxiliary table for %s\n", path_to_snp_dict);
        return 3;
    }
    
    for(i = 0; i < dict_siz; ++i) {
        packed = read_pac(fp);
        kmer_count = read_count(fp); //kmer counter
        buffer = read_flag(fp); //SNP info
        buffer = read_flag(fp); //AMBIGUOUS FLAG
        buffer = read_flag(fp); //ref freq
        buffer = read_flag(fp); //alt freq
        decode_kmer(packed, kmer);
        printf(">k\n%s\n", kmer); //print the kmer only
    }
    
    fclose(fp);
    return 0;
}

int main(int argc, char** argv)
{
    pac_t counter = 0;
    for(int i=1; i < argc; ++i)
    {
        if(strcmp(argv[i], "-r") == 0)
        {
            //reference file
            ++i;
            if(print_reference(argv[i], &counter)) return -1;
        } 
        else if(strcmp(argv[i], "-s") == 0)
        {
            //SNP file
            ++i;
            if(print_snp_dict(argv[i], &counter)) return -1;
        }
        else if(strcmp(argv[i], "-q") == 0)
        {
            //quartz dict
            ++i;
            if(print_quartz_dict(argv[i], &counter)) return -1;
        }
        else
        {
            fprintf(stderr, "Usage:\n");
            fprintf(stderr, "-r <path_to_ref_database> to print a LAVA reference database (the equivalent to the reference FASTA)\n");
            fprintf(stderr, "-s <path_to_snp_database> to print a LAVA SNP database (the database generated from the SNPs list)\n");
            fprintf(stderr, "-q <path_to_quartz_database> to print a Quartz database\n");
        }
    }
}
