#pragma once

#define ED_ASSERT(x) while(!(x)) __builtin_trap()
