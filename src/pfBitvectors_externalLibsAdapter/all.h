// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#ifdef PFBITVECTORS_USE_PASTA
#include "bitvectors/Pasta_FlatRank.h"
#include "bitvectors/Pasta_WideRank.h"
#endif

#ifdef PFBITVECTORS_USE_SDSL
#include "bitvectors/sdsl_v.h"
#include "bitvectors/sdsl_v5.h"
#endif

#ifdef PFBITVECTORS_USE_SUX
#include "bitvectors/sux_Rank9.h"
#endif

#ifdef PFBITVECTORS_USE_SDSL
#include "strings/Sdsl_wt_bldc.h"
#include "strings/Sdsl_wt_epr.h"
#endif

#ifdef PFBITVECTORS_USE_AWFMINDEX
#include "strings/AWFMIndex.h"
#endif
