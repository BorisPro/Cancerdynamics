#ifndef EVENTS_H
#define EVENTS_H

#include <vector>

enum EventType { BIRTH, DEATH, SWITCH, PRODUCTION, KILL, NONE };

class EventClass
{
public:
    EventClass();
    void clearEvents();

public:
    double eventTime;
    std::size_t chosenTrait;
    EventType type;
};

#endif // EVENTS_H
