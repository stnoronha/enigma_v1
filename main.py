 This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


from enigmatoolbox.datasets import load_summary_stats


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.




# Get case-control cortical thickness and surface area tables

#CT = sum_stats['CortThick_case_vs_controls']
#SA = sum_stats['CortSurf_case_vs_controls']

# Extract Cohen's d values

#CT_d = CT['d_icv']
#SA_d = SA['d_icv']

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    sum_stats = load_summary_stats('22q')
    CT = sum_stats['CortThick_case_vs_controls']
    SA = sum_stats['CortSurf_case_vs_controls']

    # Extract Cohen's d values

    CT_d = CT['d_icv']
    SA_d = SA['d_icv']

    print(CT_d)
