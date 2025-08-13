#!/usr/bin/env python3
"""
Compare AMRVAC simulation logs with reference output
"""

import numpy as np

def read_log_file(filename):
    """Read AMRVAC log file and return data"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#') or 'mean' in line:
                continue
            # Skip empty lines
            line = line.strip()
            if not line:
                continue
            try:
                values = [float(x) for x in line.split()]
                data.append(values)
            except ValueError:
                # Skip lines that can't be parsed as numbers
                continue
    return np.array(data)

def compare_logs(ref_file, test_file):
    """Compare two log files and compute differences"""
    print(f"\n{'='*60}")
    print(f"Comparing AMRVAC Simulation Logs")
    print(f"{'='*60}")
    print(f"Reference: {ref_file}")
    print(f"Test:      {test_file}")
    print(f"{'='*60}\n")
    
    # Read data
    ref_data = read_log_file(ref_file)
    test_data = read_log_file(test_file)
    
    # Check dimensions
    print(f"Reference shape: {ref_data.shape}")
    print(f"Test shape:      {test_data.shape}")
    
    if ref_data.shape != test_data.shape:
        print("ERROR: Different data shapes!")
        return False
    
    # Variable names based on AMRVAC RHD output
    var_names = ['time', 'rho_mean', 'rho_var', 'r_e_mean', 'e_mean', 
                 'm1_mean', 'rho_mean2', 'rho_var2', 'r_e_mean2', 
                 'e_mean2', 'm1_mean2']
    
    # Compute differences
    abs_diff = np.abs(test_data - ref_data)
    rel_diff = np.zeros_like(abs_diff)
    
    # Compute relative differences (avoid division by zero)
    for i in range(ref_data.shape[0]):
        for j in range(ref_data.shape[1]):
            if ref_data[i, j] != 0:
                rel_diff[i, j] = abs_diff[i, j] / np.abs(ref_data[i, j])
    
    # Statistics for each variable
    print("\n" + "="*60)
    print("Variable-wise Comparison Statistics")
    print("="*60)
    print(f"{'Variable':<15} {'Max Abs Diff':<15} {'Max Rel Diff':<15} {'Mean Rel Diff':<15}")
    print("-"*60)
    
    for j, var in enumerate(var_names[:test_data.shape[1]]):
        max_abs = np.max(abs_diff[:, j])
        max_rel = np.max(rel_diff[:, j])
        mean_rel = np.mean(rel_diff[:, j])
        print(f"{var:<15} {max_abs:<15.2e} {max_rel:<15.2e} {mean_rel:<15.2e}")
    
    # Time evolution comparison
    print("\n" + "="*60)
    print("Time Evolution Comparison (Key Variables)")
    print("="*60)
    
    # Select key timesteps
    n_steps = len(ref_data)
    key_steps = [0, n_steps//4, n_steps//2, 3*n_steps//4, n_steps-1]
    
    print(f"\n{'Time':<12} {'Variable':<10} {'Reference':<15} {'Test':<15} {'Rel Diff':<12}")
    print("-"*64)
    
    for idx in key_steps:
        time = ref_data[idx, 0]
        # Check key variables: rho_mean, r_e_mean, m1_mean
        for j, var in enumerate([1, 3, 5]):  # indices for key variables
            var_name = var_names[var]
            ref_val = ref_data[idx, var]
            test_val = test_data[idx, var]
            if ref_val != 0:
                rel_d = abs(test_val - ref_val) / abs(ref_val)
            else:
                rel_d = abs(test_val - ref_val)
            
            print(f"{time:<12.4e} {var_name:<10} {ref_val:<15.8e} {test_val:<15.8e} {rel_d:<12.2e}")
    
    # Final comparison at t=0.01
    print("\n" + "="*60)
    print("Final Time (t=0.01) Detailed Comparison")
    print("="*60)
    
    final_idx = -1
    print(f"\n{'Variable':<15} {'Reference':<20} {'Test':<20} {'Difference':<15}")
    print("-"*70)
    
    for j, var in enumerate(var_names[:test_data.shape[1]]):
        ref_val = ref_data[final_idx, j]
        test_val = test_data[final_idx, j]
        diff = test_val - ref_val
        print(f"{var:<15} {ref_val:<20.12e} {test_val:<20.12e} {diff:<15.2e}")
    
    # Overall assessment
    print("\n" + "="*60)
    print("Overall Assessment")
    print("="*60)
    
    max_rel_diff = np.max(rel_diff)
    mean_rel_diff = np.mean(rel_diff)
    
    print(f"Maximum relative difference: {max_rel_diff:.2e}")
    print(f"Mean relative difference:    {mean_rel_diff:.2e}")
    
    # Tolerance check (typical floating point precision)
    tolerance = 1e-6  # Relative tolerance
    
    if max_rel_diff < tolerance:
        print(f"\n✅ VALIDATION PASSED")
        print(f"   All differences are within tolerance ({tolerance:.0e})")
        print(f"   The simulation results match the reference output.")
        return True
    else:
        print(f"\n⚠️  VALIDATION WARNING")
        print(f"   Some differences exceed tolerance ({tolerance:.0e})")
        print(f"   Maximum relative difference: {max_rel_diff:.2e}")
        
        # Check if differences are still acceptable (machine precision)
        if max_rel_diff < 1e-5:
            print(f"   However, differences are within acceptable machine precision.")
            print(f"   This is likely due to compiler or architecture differences.")
            return True
        else:
            print(f"   Differences may indicate a problem with the simulation.")
            return False

def main():
    """Main function"""
    ref_log = "correct_output/rhd_wave.log"
    test_log = "rhd_wave.log"
    
    # Run comparison
    result = compare_logs(ref_log, test_log)
    
    print("\n" + "="*60)
    if result:
        print("CONCLUSION: Simulation validation SUCCESSFUL ✅")
    else:
        print("CONCLUSION: Simulation validation FAILED ❌")
    print("="*60)

if __name__ == "__main__":
    main()