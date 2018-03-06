#include "BaseRecord.hpp"

namespace Tracing{

BaseRecord::BaseRecord(std::string ID) noexcept
{
    this->id = ID;
    ctor();
}

std::string BaseRecord::getId() const
{
    return id;
}

std::string BaseRecord::getName() const
{
    return name;
}

std::string BaseRecord::getDescription() const
{
    return description;
}

void BaseRecord::setName(std::string nname)
{
    this->name = nname;
}

void BaseRecord::setDescription(std::string ndescription)
{
    this->description = ndescription;
}

void BaseRecord::ctor()
{
    name = "";
    description = "";
}

BaseRecord::~BaseRecord()
{
    //dtor
}

}//Tracing
