import matplotlib.pyplot as plt

# Data for the bar diagrams
meshes = ['bunny', 'eba', 'crank_pin']

# all data is million queries in milliseconds

# on apple m2
mantis_times = [149.962875, 616.285834, 542.318333]
p2m_times = [320.562167, 1389.114166, 812.242750]
fcpw_times = [645.630916, 2835.470958, 1641.606958]

# on windows with MinGW and Ryzen 9 7950X cpu
#mantis_times = [93.728900, 330.634800, 335.202500]
#p2m_times = [400.812200, 1533.474400, 1041.552700]
#fcpw_times = [1684.176200, 7286.450900, 4191.863200]

# Convert milliseconds to microseconds for 1 query
mantis_times_us = [time / 1000 for time in mantis_times]
p2m_times_us = [time / 1000 for time in p2m_times]
fcpw_times_us = [time / 1000 for time in fcpw_times]

# Calculate the number of times slower than mantis
p2m_slower = [p2m / mantis for p2m, mantis in zip(p2m_times_us, mantis_times_us)]
fcpw_slower = [fcpw / mantis for fcpw, mantis in zip(fcpw_times_us, mantis_times_us)]

# Define new color scheme with switched colors for Mantis and P2M
colors = ['lightgreen', 'skyblue', 'salmon']

# Create bar diagrams
fig, ax = plt.subplots(figsize=(10, 6))
x = range(len(meshes))

# Plotting the bars with new color scheme
bar_width = 0.2
ax.bar([p - bar_width for p in x], mantis_times_us, width=bar_width, label='Mantis', color=colors[0])
ax.bar(x, p2m_times_us, width=bar_width, label='P2M', color=colors[1])
ax.bar([p + bar_width for p in x], fcpw_times_us, width=bar_width, label='FCPW', color=colors[2])

# Adding the slowdown factors on top of the bars for P2M and FCPW

# Adding the slowdown factors on top of the bars for P2M and FCPW with adjusted settings
for i in x:
    ax.text(i - bar_width, mantis_times_us[i] + max(fcpw_times_us) * 0.0, '1x', ha='center', va='bottom')
    ax.text(i, p2m_times_us[i] + max(fcpw_times_us) * 0.0, f'{p2m_slower[i]:.1f}x', ha='center', va='bottom')
    ax.text(i + bar_width, fcpw_times_us[i] + max(fcpw_times_us) * 0.0, f'{fcpw_slower[i]:.1f}x', ha='center', va='bottom')

# Setting the x-axis labels
ax.set_xticks(x)
ax.set_xticklabels(meshes)

# Setting the y-axis label
ax.set_ylabel('time in μs')

# Adding a legend
ax.legend()

ax.grid(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Show the updated plot
plt.show()