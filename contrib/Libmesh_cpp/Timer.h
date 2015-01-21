#ifndef _TIMER_H
#define _TIMER_H

#include <time.h>
#include <sys/time.h>

class Timer {

  struct timeval _t;
  double _t0, _t1, _spent;
  bool _paused;

public:

  Timer() : _spent(0.0), _paused(true) {
  }

  void start() {
    gettimeofday(&_t, NULL);
    _t0 = (double)_t.tv_sec + (double)_t.tv_usec * 1e-6;
    _spent = 0.0;
    _paused = false;
  }

  double elapsed() {
    if (_paused) return _spent;
    
    gettimeofday(&_t, NULL);
    _t1 = (double)_t.tv_sec + (double)_t.tv_usec * 1e-6;
    return _spent + (_t1 - _t0);
  }

  void pause() {
    if (_paused) return;
    _spent = elapsed();
    _paused = true;
  }

  void resume() {
    if (!_paused) return;

    gettimeofday(&_t, NULL);
    _t0 = (double)_t.tv_sec + (double)_t.tv_usec * 1e-6;
    _paused = false;
  }
  
};

#endif
