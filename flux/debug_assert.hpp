// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#ifndef FLUX_DEBUG_ASSERT_HPP_
#define FLUX_DEBUG_ASSERT_HPP_

#include <mpi.h>
#include <iostream>

#ifdef NDEBUG
#define DEBUF_ASSERT_DISABLE
#define FLUX_DEBUG_ASSERT_LEVEL 1
#endif

#ifndef FLUX_DEBUG_ASSERT_LEVEL
#define FLUX_DEBUG_ASSERT_LEVEL 9
#endif

// This header contains the foonathan debug_assert code that we will be using.
// Need to be after other defines to perform preprocessor removal of
// DEBUG_ASSERT
#include <debug_assert/debug_assert.hpp>

  struct global_assert : debug_assert::set_level<FLUX_DEBUG_ASSERT_LEVEL> {
  static void handle(const debug_assert::source_location& loc,
                     const char* expression,
                     std::string message = "") noexcept {
    std::cerr << "Assertion failure '" << loc.file_name << ':'
              << loc.line_number << ":\n";
    std::cerr << "Failing expression: " << expression;
    if (message != "") {
      std::cerr << "\nMessage: " << message;
    }
    std::cerr << '\n';
    int flag;
    MPI_Initialized(&flag);
    if (flag == 1) {
      // TODO, make this not MPI_COMM_WORLD
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  };

    namespace DebugLevel {
// Largest amount of debugging. Used for full debug. Can include expensive
// checks.
using FULL = debug_assert::level<9>;

// Cheap debugging checks, such as checking bounds on an object or index. A
// check that is not expected to greatly increase computational cost. Allows
// debug runs to be done more efficiently if full option not needed.
using CHEAP = debug_assert::level<5>;

// To never be turned off. Used as critical or fatal failures.
using ALWAYS = debug_assert::level<1>;
    }  // namespace DebugLevel


// Custom macro for fatal errors
// Reuse some DEBUG_ASSERT magic
#define FLUX_FATAL_ERROR(Expr, ...)                                     \
  static_cast<void>(debug_assert::detail::do_assert(                     \
      [&]() noexcept { return Expr; }, DEBUG_ASSERT_CUR_SOURCE_LOCATION, \
      #Expr, global_assert{}, debug_assert::level<1>{}, __VA_ARGS__))

#endif  // FLUX_DEBUG_ASSERT_HPP_
