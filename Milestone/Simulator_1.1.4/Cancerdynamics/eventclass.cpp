#include "eventclass.h"

EventClass::EventClass()
{
}

void EventClass::clearEvents()
{
    eventTime = 0;
    chosenTrait = 0;
    type = NONE;
//    isBirth = false;
}
