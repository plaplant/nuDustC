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

#include <plog/Log.h>

#include <plog/Initializers/RollingFileInitializer.h>

#include <plog/Formatters/TxtFormatter.h>

#include <plog/Appenders/ColorConsoleAppender.h>

//#include <plog/Appenders/DynamicAppender.h>

enum

{

  DetailLog_Root = 1,

  DetailLog_Rank = 2

};

static plog::ColorConsoleAppender<plog::TxtFormatter> rootAppender;

static plog::ColorConsoleAppender<plog::TxtFormatter> rankAppender;

#define SCLOG_0(rank) if (rank == 0) PLOGI_(DetailLog_Root)

#define SCLOG_A(rank) PLOGI_(DetailLog_Root) << "[" << rank<< "] "

