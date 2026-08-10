import numpy as np
import scipy.stats as stats

# 1. Input your summary statistics for the 10 samples
SPd = np.array([
    0.410848, 0.356645, 0.528239, 0.477902, 0.448119,
    0.432432, 0.357712, 0.600614, 0.485576, 0.515038
])

stdevs = np.array([
    0.031534, 0.04823, 0.047206, 0.058493, 0.029654, 
    0.03321, 0.048378, 0.053729, 0.059462, 0.034133
])

# The number of observations that made up each mean 
ns = np.array([3] * 10) 

names = [
    "R1a", "R2a", "R3a", "R4a", "R5a", 
    "R1b", "R2b", "R3b", "R4b", "R5b"
]

# Set significance level
alpha = 0.05

for i in range(10):
    # Wrap the slices in an extra set of parentheses so they form a tuple
    means = np.concatenate((SPd[:i], SPd[i+1:]))
    stds =  np.concatenate((stdevs[:i], stdevs[i+1:]))
    n = np.concatenate((ns[:i], ns[i+1:]))

    # 2. Calculate the Grand Mean
    grand_mean = np.sum(n * means) / np.sum(n)

    # 3. Calculate Sum of Squares Between (SSB)
    ssb = np.sum(n * (means - grand_mean)**2)
    df_between = len(means) - 1

    # 4. Calculate Sum of Squares Within (SSW)
    ssw = np.sum((n - 1) * stds**2)
    df_within = np.sum(n) - len(means)

    # 5. Calculate Mean Squares
    msb = ssb / df_between
    msw = ssw / df_within

    # 6. Calculate the F-statistic and P-value
    f_stat = msb / msw
    p_value = 1 - stats.f.cdf(f_stat, df_between, df_within)

    # 7. Evaluate significance
    is_significant = "Significant" if p_value < alpha else "Not significant"

    print(f"Without {names[i]:<4} & P-Value: {p_value:.4e} & Result: {is_significant}")