/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#pragma once

#include<string>
#include<chrono>
#include <unistd.h>

using tpoint = std::chrono::time_point<std::chrono::steady_clock>;
using trec = std::map<std::string, std::vector<double>>;

inline trec benchmark_map;

struct Timer
{
    std::string fn;
    tpoint start;
    Timer(std::string fn)
        : fn(std::move(fn)), start(std::chrono::steady_clock::now())
    {
    }
    ~Timer()
    {
        const auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//       printf("%s: function=%s; elasepd=%f ms\n", title.c_str(), fn.c_str(), elapsed / 1000.0);
        benchmark_map[fn].push_back(static_cast<double>(elapsed));
    }
  
};


#ifndef ENABLE_BENCHMARK
static constexpr inline void dummy_fn() { }
#define START_BENCHMARK_TIMER(...) dummy_fn()
#else
#define START_BENCHMARK_TIMER(title) Timer timer(title)
#endif

template<class F, typename ...Args>
auto time_fn(std::string fn_name, F&& fn, Args&&... args ) {
  START_BENCHMARK_TIMER(fn_name);
  return fn(std::forward<Args>(args)...);
}
