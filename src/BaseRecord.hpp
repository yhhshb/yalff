#ifndef BASERECORD_HPP
#define BASERECORD_HPP

#include <string>
#include <vector>

namespace Tracing{

class BaseRecord
{
    public:
        BaseRecord(const BaseRecord&) = default;

        BaseRecord(BaseRecord&&) = default;

        BaseRecord(std::string id) noexcept;

        virtual std::string getId() const;

        virtual std::string getName() const;
        virtual std::string getDescription() const;

        virtual void setName(std::string name);
        virtual void setDescription(std::string description);

        virtual ~BaseRecord();

    protected:
        void ctor();

        std::vector<std::string> dbxrefs;
        //std::vector<SeqFeature> features;
        std::vector<std::string> annotations;

        std::string id;
        std::string name;
        std::string description;
};

}//Tracing

#endif // BASERECORD_HPP
