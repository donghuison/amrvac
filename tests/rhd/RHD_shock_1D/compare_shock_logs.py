#!/usr/bin/env python3
"""
Compare RHD shock simulation log with reference output
"""

import re

def extract_final_values(log_file):
    """Extract final values from regression test log file"""
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    # This is a regression test log with time series data
    # Format: time mean(w) mean(w**2) for each variable
    # Find the last data line
    last_data_line = None
    for line in reversed(lines):
        line = line.strip()
        if line and not line.startswith('#'):
            # Check if it's a data line (starts with a number)
            try:
                float(line.split()[0])
                last_data_line = line
                break
            except (ValueError, IndexError):
                continue
    
    if last_data_line:
        # Parse the line
        parts = last_data_line.split()
        if len(parts) >= 11:  # Has all columns
            return {
                'time': float(parts[0]),
                'mean_rho': float(parts[1]),
                'mean_mom1': float(parts[2]),
                'mean_e': float(parts[3]),
                'mean_r_e': float(parts[4]),
                'mean2_rho': float(parts[5]),
                'mean2_mom1': float(parts[6]),
                'mean2_e': float(parts[7]),
                'mean2_r_e': float(parts[8])
            }
    return None

def compare_logs(test_log, ref_log):
    """Compare test log with reference log"""
    print("="*60)
    print("RHD Shock Simulation Log Comparison")
    print("="*60)
    
    test_values = extract_final_values(test_log)
    ref_values = extract_final_values(ref_log)
    
    if not test_values:
        print(f"Could not extract values from {test_log}")
        return
    
    if not ref_values:
        print(f"Could not extract values from {ref_log}")
        return
    
    print(f"\nTest log: {test_log}")
    print(f"Reference log: {ref_log}")
    print("\n" + "-"*40)
    
    # Compare final time
    print(f"\nFinal time:")
    print(f"  Test:      {test_values['time']:.6E}")
    print(f"  Reference: {ref_values['time']:.6E}")
    diff_time = abs(test_values['time'] - ref_values['time'])
    print(f"  Difference: {diff_time:.6E}")
    
    # Compare mean values at final time
    print(f"\nMean values at t={test_values['time']:.2f}:")
    
    variables = [
        ('Density (rho)', 'mean_rho'),
        ('Momentum (mom1)', 'mean_mom1'),
        ('Energy (e)', 'mean_e'),
        ('Radiation Energy (r_e)', 'mean_r_e')
    ]
    
    max_rel_diff = 0.0
    for var_name, var_key in variables:
        test_val = test_values[var_key]
        ref_val = ref_values[var_key]
        diff = abs(test_val - ref_val)
        rel_diff = diff / abs(ref_val) * 100 if ref_val != 0 else 0
        max_rel_diff = max(max_rel_diff, rel_diff)
        
        print(f"\n  {var_name}:")
        print(f"    Test:      {test_val:.6E}")
        print(f"    Reference: {ref_val:.6E}")
        print(f"    Rel. Diff: {rel_diff:.6E}%")
    
    # Summary
    print("\n" + "="*60)
    if max_rel_diff < 1e-8:
        print("✅ VALIDATION PASSED: Exact match with reference!")
    elif max_rel_diff < 1e-4:
        print("✅ VALIDATION PASSED: Results match within numerical precision")
    elif max_rel_diff < 0.1:
        print("✅ VALIDATION PASSED: Results match within 0.1% tolerance")
    else:
        print(f"⚠️  Maximum relative difference: {max_rel_diff:.2E}%")
    print("="*60)

def main():
    """Main function to run the log comparisons"""
    # Compare with reference
    test_log = 'rhd_shock.log'
    ref_log = 'correct_output/rhd_shock.log'
    
    compare_logs(test_log, ref_log)
    
    # Also compare with IMEX_SP_tvdlf_ko variant if it exists
    print("\n")
    ref_log2 = 'correct_output/rhd_shock_IMEX_SP_tvdlf_ko.log'
    print("Comparing with IMEX_SP_tvdlf_ko variant:")
    compare_logs(test_log, ref_log2)

if __name__ == "__main__":
    main()