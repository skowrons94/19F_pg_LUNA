/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorManager.hh"
#include "RunAction.hh"
#include "EventAction.hh"

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorManager());
  SetUserAction(new RunAction());
  SetUserAction(new EventAction());
}
