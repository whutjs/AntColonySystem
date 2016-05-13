#ifndef TIMER_H
#define TIMER_H
typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;

void start_timers(void);
double elapsed_time(TIMER_TYPE type);

#endif