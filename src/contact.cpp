#include "contact.h"

namespace mc_impact
{

McContact::McContact(const ContactParams & inputParams) : contactParams_(inputParams)
{

  updateCWC();
}

void McContact::updateCWC()
{
  double X = getContactParams().halfX;
  double Y = getContactParams().halfY;
  double mu = getContactParams().frictionCoe;

  CWC_.resize(16, 6);
  CWC_.setZero();
  CWC_ <<
      // mx,  my,  mz,  fx,  fy,            fz,
      0,
      0, 0, -1, 0, -mu, 0, 0, 0, +1, 0, -mu, 0, 0, 0, 0, -1, -mu, 0, 0, 0, 0, +1, -mu, -1, 0, 0, 0, 0, -Y, +1, 0, 0, 0,
      0, -Y, 0, -1, 0, 0, 0, -X, 0, +1, 0, 0, 0, -X, +mu, +mu, -1, -Y, -X, -(X + Y) * mu, +mu, -mu, -1, -Y, +X,
      -(X + Y) * mu, -mu, +mu, -1, +Y, -X, -(X + Y) * mu, -mu, -mu, -1, +Y, +X, -(X + Y) * mu, +mu, +mu, +1, +Y, +X,
      -(X + Y) * mu, +mu, -mu, +1, +Y, -X, -(X + Y) * mu, -mu, +mu, +1, -Y, +X, -(X + Y) * mu, -mu, -mu, +1, -Y, -X,
      -(X + Y) * mu;
}

const McContact & McContactSet::getContact(const std::string & name)
{
  auto opt = contacts_.find(name);
  if(opt != (contacts_.end()))
  {
    return opt->second;
  }
  else
  {
    throw std::runtime_error(std::string("McContactSet::getContact: '") + name + std::string("' is not found."));
  }
}

bool McContactSet::addContact(const ContactParams & inputParams)
{
  auto opt = contacts_.find(inputParams.surfaceName);

  if(opt != (contacts_.end()))
  {
    std::cout << "Endeffector " << inputParams.surfaceName << " already exists." << std::endl;
    return false;
  }

  contacts_.insert({"test", {inputParams}});

  std::cout << "McContactSet: Adding end-effector contact: " << inputParams.surfaceName << std::endl;

  return true;
}

} // namespace mc_impact
