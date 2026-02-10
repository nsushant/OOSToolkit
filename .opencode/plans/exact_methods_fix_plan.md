# Exact Methods Fix Implementation Plan

## Issues Identified (Based on Direct Dynamic Programming Call)

After analyzing `sim.cpp` and how you actually call the dynamic programming function (line 182):
```cpp
int result = dynamic_program_fixed_tasksize_Tfixed(schedule_base, tstep_size, 0.0, schedule_base.blocks.size(), Tmax, simfile,lookup);
```

The issues are:

1. **Constraint validation incorrectly applied to Block 0** - causing the error message
2. **Schedule modifications not persisted through recursive calls** - fragmented state management across recursion levels

## Root Cause Analysis

### Problem 1: Constraint Validation (Line 197)
- Block 0 (starting depot) should not have `departure_time <= arrival_constraint` validation
- This causes: "Constraint violation in decrement_optimization:Block 0 departure_time (6000.00) > arrival_constraint (0.00)"
- Block 0 typically has `arrival_constraint = 0.0` and should be excluded from this check

### Problem 2: Schedule Persistence in Dynamic Programming
The `dynamic_program_fixed_tasksize_Tfixed` function has structural issues:

1. **Fragmented State Management**: Each recursion level creates its own `ScheduleStateManager` instance
2. **No Final State Consolidation**: Changes committed locally aren't properly propagated back
3. **Recursion Direction Mismatch**: `b_reached` increments forward, but optimization runs in reverse

## Specific Code Changes Required

### 1. Fix Constraint Validation (Line 197)
**File**: `src/Exact_methods.cpp`
**Function**: `validate_schedule_constraints`
**Change**: Skip arrival constraint check for Block 0

```cpp
// Current (line 197):
if (block.departure_time >= block.arrival_constraint ) {

// Fixed:
if (i > 0 && block.departure_time >= block.arrival_constraint ) {
```

### 2. Fix Schedule State Management (Lines 603, 638)
**File**: `src/Exact_methods.cpp`
**Function**: `dynamic_program_fixed_tasksize_Tfixed`

**Problem**: Each recursion creates separate `ScheduleStateManager` instances, causing fragmented state

**Solution A**: Use single ScheduleStateManager across all recursion levels
**Solution B**: Ensure final state commitment after recursive calls complete

**Option A (Preferred) - Refactor to use shared state manager:**
- Move ScheduleStateManager creation outside the recursive function
- Pass it as a parameter to recursive calls
- Ensure single point of state management

**Option B - Ensure final commitment:**
- Add final state consolidation after recursive calls
- Ensure the schedule reference maintains all committed changes

## Implementation Order

1. **Fix constraint validation** (quick fix - eliminates the error)
2. **Fix schedule state management** (main fix - ensures modifications persist)

## Expected Results

After these changes:
- ✅ No more Block 0 constraint violation errors
- ✅ Schedule modifications from dynamic programming will be persisted to caller's `schedule_base`
- ✅ The `view_schedule(schedule_base)` call will show the optimized schedule
- ✅ `deltavtotcalc(schedule_base)` will calculate the optimized deltaV

## Testing Strategy

1. **Immediate**: Fix constraint validation and verify error disappears
2. **Verify**: Check that `schedule_base` is modified after the dynamic programming call
3. **Confirm**: Compare schedule before/after to ensure improvements are applied
4. **Integration**: Run with your scaling tests to ensure consistent results

**Note**: The `branch` and `branch_and_bound` functions are not relevant since you're calling the dynamic programming function directly.