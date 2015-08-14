"""
Functions and such to parse and plot timing output from SELFE.

jlopez 2013-05-18
"""
#-------------------------------------------------------------------------------
# Import
#-------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import glob
import os

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------
TIMER_IDS = {0 : 'Total',
             1 : 'Initialize',
             2 : 'Timestep',
             3 : 'Forcing & Prep', 
             4 : 'Backtrack',
             5 : 'Turb Closure',
             6 : 'Matrix Prep',
             7 : 'Wave-Continuity',
             8 : 'Momentum',
             9 : 'Transport', 
             10 : 'Levels',
             11 : 'Cons Check',
             12 : 'Output',
             13 : 'Hotstart',}

TIMER_IDS_L = ['Total',
              'Initialize',
              'Timestep',
              'Forcing & Prep', 
              'Backtrack',
              'Turb Closure',
              'Matrix Prep',
              'Wave-Continuity',
              'Momentum',
              'Transport', 
              'Levels',
              'Cons Check',
              'Output',
              'Hotstart',]

COL_PLOT = ['Forcing & Prep',
            'Backtrack',
            'Turb Closure',
            'Matrix Prep',
            'Wave-Continuity',
            'Momentum',
            'Transport', 
            'Levels',
            'Cons Check',
            'Output',
            'Hotstart',]

DESC_STATS_COLS = { 0 : 'ID',
                    1 : 'Comp',
                    2 : 'Min comp',
                    3 : 'Min comp rank',
                    4 : 'Max comp',
                    5 : 'Max comp rank',
                    6 : 'Comm',
                    7 : 'Min comm',
                    8 : 'Min comm rank',
                    9 : 'Max comm',
                    10 : 'Max comm rank', }

DESC_STATS_COLS_L = ['ID',
                     'Comp',
                     'Min comp ',
                     'Min comp rank',
                     'Max comp',
                     'Max comp rank',
                     'Comm',
                     'Min comm',
                     'Min comm rank',
                     'Max comm',
                     'Max comm rank', ]
#-------------------------------------------------------------------------------
# Functions to read file 
#-------------------------------------------------------------------------------
def read_file(path):
  """ Function to read a SELFE timer.out file 
  
  Parameters:
  -----------
    * path - String of path to file 

  Returns:
  --------
    * times - Dict of data structure with times 
  """
  times = {}
  with open(path) as f:
    _read_header(f)
    times['desc_stats'] = _read_desc_stats(f)
    times['comp_times'] = _read_times_section(f) 
    times['comm_times'] = _read_times_section(f)

  return times

def _read_header(file):
  """ Reads the header for the file. """
  # Header + 14 categories + 2 blank lines
  for i in range(1 + 14 + 2):
    file.readline()


def _read_desc_stats(file):
  """ Reads the descriptive statistics section. """
  # 2 header lines + 14 categories + 2 blank lines
  desc = []
  for i in range(2): 
    file.readline()
  for i in range(14): 
    tmp = file.readline().split() 
    desc.append(dict((DESC_STATS_COLS[j], float(v)) for (j, v) in enumerate(tmp)))   
  for i in range(2):
    file.readline()

  desc_df = pd.DataFrame(desc, index=TIMER_IDS_L) 
  desc_df.index.name = 'Sections'

  return desc_df

def _read_times_section(file):
  """ Reads the computation or times section. """
  rank = -1
  rank_times = [] 
   # 3 header lines + rank info + 2 blank lines
  for i in range(3):
    file.readline()
  while True:
    tmp = file.readline().split() 
    if not tmp: 
      break
    rank += 1
    # [2:] to skip useless # and rank from lines
    rank_times.append(dict((TIMER_IDS[i], float(v)) for (i, v) in enumerate(tmp[2:])))
  for i in range(1):
      file.readline()

  rank_times_df = pd.DataFrame(rank_times, columns=TIMER_IDS_L)
  rank_times_df.index.name = 'Rank'
  rank_times_df.columns.name = 'Section'

  return rank_times_df


#-------------------------------------------------------------------------------
# Functions for working with ensembles 
#-------------------------------------------------------------------------------
def load_timer_set(dir, file_proto='timer*'):
  """Loads a set of timer.out files and returns list of timers
  
  Parameters
  ----------
    * dir - Directory where timers are stored
    * file_proto - String to match timer file names against

  Returns
  -------
    * timers - List of timers
  """
  timers = []
  for f in glob.glob(file_proto):
    fp = os.path.join(dir, f)
    timers.append(read_file(fp))

  return timers

def order_timers_by_proc(timers):
  """Return timers by number of ranks. Useful for scaling tests.

  Parameters
  ----------
    * timers - List of timers

  Returns
  -------
    * ranks - List of ranks used for test
    * sorted - List of timers ordered 
  """
  ranks = []
  for t in timers:
    ranks.append(len(t['comp_times']['Total']))
  rank2timer_map = dict((ranks[o], o) for o in range(len(ranks)))
  ranks.sort() 

  sorted = []
  for r in ranks:
    sorted.append(timers[rank2timer_map[r]])     

  return ranks, sorted
#-------------------------------------------------------------------------------
# Functions for calculating statistics  
#-------------------------------------------------------------------------------
def calc_section_percents(times, time_type):
  """Calculates the percentage of time spent per section per rank.
   
  Parameters: 
  -----------
  * times - Dict form read_times()
  * time_type - Computation ('comp') or communicaiton ('comm')

  Returns:
  --------
  * percents - DataFrame of percent of time by rank and section

  """
  if time_type == 'comp':
    key = 'comp_times'
  elif time_type == 'comm':
    key = 'comm_times'
  else:
    raise Exception('Type not supported: '+time_type)

  percents = times[time_type].div(times[time_type].sum(1), axis=0)

  return percents

def calc_total_time(times):
  """Calculate the total time for a run."""
  total = times['comp_times']['Total'].max() + times['comm_times']['Total'].max()

  return total


#-------------------------------------------------------------------------------
# Functions to plot data 
#-------------------------------------------------------------------------------
def plot_times(times, plot_type, path):
  """Creates a stacked barplot of computation times by rank.

  Parameters:
  -----------
  * times - Dict output from read_timer()
  * plot_type - Computation ('comp') or communication ('comm')
  * path - Path/name of output image
  """
  if plot_type == 'comp':
    key = 'comp_times'
    title = 'Computation Times'
  elif plot_type == 'comm':
    key = 'comm_times'
    title = 'Communication Times'
  else:
    raise Exception('Type not supported: '+plot_type)

  sections = times[key][COL_PLOT]
  ax = sections.plot(kind='barh', stacked=True, legend=False, title=title)
  patches, labels = ax.get_legend_handles_labels()
  ax.set_position([0.1, 0.1, 0.5, 0.8])
  ax.set_xlabel('Time (s)')
  ax.legend(patches, labels, loc='center left', bbox_to_anchor=(1.0, 0.5))
  figure = plt.gcf()
  figure.savefig(path)


def plot_all_times(times, path):
  """Creates a plot of computation and communication times by rank.

  Parameters:
  -----------
  * times - Dict output from read_timer()
  * path - Path/name of output image
  """
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,10))
  sections1 = times['comp_times'][COL_PLOT]
  sections2 = times['comm_times'][COL_PLOT]
  ax1 = sections1.plot(kind='bar', legend=False, ax=axes[0]) 
  ax2 = sections2.plot(kind='bar', legend=COL_PLOT, ax=axes[1])
  axes[0].set_title('Computation Time')
  axes[1].set_title('Communication Time')
  ax1.set_ylabel('Time (seconds)')
  ax2.set_ylim(ax1.get_ylim())
  ax2.set_yticklabels('')

  plt.tight_layout()
  fig.savefig(path)


def plot_all_times_stacked(times, path):
  """Creates a stacked barplot of computation and communication times by rank.

  Parameters:
  -----------
  * times - Dict output from read_timer()
  * path - Path/name of output image
  """
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,10))
  sections1 = times['comp_times'][COL_PLOT]
  sections2 = times['comm_times'][COL_PLOT]
  ax1 = sections1.plot(kind='barh', stacked=True, legend=False, ax=axes[0]) 
  ax2 = sections2.plot(kind='barh', stacked=True, legend=COL_PLOT, ax=axes[1])
  axes[0].set_title('Computation Time')
  axes[1].set_title('Communication Time')
  ax1.set_xlabel('Time (seconds)')
  ax2.set_xlabel('Time (seconds)')

  plt.tight_layout()
  fig.savefig(path)


def plot_normalized_times_stacked(times, path):
  """Creates a normalized stacked barplot of computation and communication 
  times by rank.

  Parameters:
  -----------
  * times - Dict output from read_timer()
  * path - Path/name of output image
  """
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,10))
  sections1 = times['comp_times'][COL_PLOT]
  sections2 = times['comm_times'][COL_PLOT]
  percents1 = sections1.div(sections1.sum(1), axis=0)
  percents2 = sections2.div(sections2.sum(1), axis=0)
  ax1 = percents1.plot(kind='bar', stacked=True, legend=COL_PLOT, ax=axes[0]) 
  ax2 = percents2.plot(kind='bar', stacked=True, legend=False, ax=axes[1])
  axes[0].set_title('Computation Time')
  axes[1].set_title('Communication Time')
  #ax1.set_xlabel('Time (seconds)')
  #ax2.set_xlabel('Time (seconds)')

  plt.tight_layout()
  fig.savefig(path)


def plot_times_averaged_by_section(times, path):
  """Creates stack barplot of computation and communication times averaged 
  across all ranks.

  Parameters:
  -----------
  * times - Dict output from read_timer()
  * path - Path/name of output image
  """
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,10))
  sections1 = times['comp_times'][COL_PLOT]
  sections2 = times['comm_times'][COL_PLOT]
  average1 = sections1.mean(0)
  average2 = sections2.mean(0)
  ax1 = average1.plot(kind='bar', stacked=True, legend=False, ax=axes[0]) 
  ax2 = average2.plot(kind='bar', stacked=True, legend=False, ax=axes[1])
  axes[0].set_title('Computation Time')
  axes[1].set_title('Communication Time')
  ax1.set_ylabel('Time (seconds)')
  #ax2.set_xlabel('Time (seconds)')

  plt.tight_layout()
  fig.savefig(path)


def plot_scaling(times_set, path, title=None, label='Test'):
  """Creates a scaling plot of ranks vs. speedup.

  You better have a rank of 1 in there or ideal isn't correct.

  Parameters:
  -----------
  * times_set - A set of times probably from load_timer_set
  * path - Path/name of output image
  """
  (r,s) = order_timers_by_proc(times_set)
  total_times = []
  for t in s:
    total_times.append(calc_total_time(t))
  ideal_times = total_times[0]/r  

  font = {'family' : 'normal',
          'size'   : 22}
  matplotlib.rc('font', **font)

  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
  ax.plot(r, ideal_times, '--b', linewidth=3, label='Ideal')
  ax.plot(r, total_times, '-r', marker='o', markersize=15, linewidth=3, label=label)

  ax.set_title(title)
  ax.set_ylabel('Time (seconds)')
  ax.set_xlabel('Processes')
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xticks(r)
  ax.set_xticklabels(r)
  plt.minorticks_off()
  ax.legend(loc='upper right', shadow=True)
  ax.grid(True)

  plt.tight_layout()
  fig.savefig(path)

def compare_averaged_by_section(times1, times2, path, name1='old', name2='new'):
  """Creates stack barplot of computation and communication times averaged 
  across all ranks.

  Parameters:
  -----------
  * times1 - Dict output from read_timer()
  * times2 - Dict output from read_times()
  * path - Path/name of output image
  """
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,10))

  # Calc average times first 
  sections1 = times1['comp_times'][COL_PLOT]
  sections2 = times1['comm_times'][COL_PLOT]
  t1_average1 = sections1.mean(0)
  t1_average2 = sections2.mean(0)

  sections1 = times2['comp_times'][COL_PLOT]
  sections2 = times2['comm_times'][COL_PLOT]
  t2_average1 = sections1.mean(0)
  t2_average2 = sections2.mean(0)

  # Plot
  avg1 = pd.DataFrame({name1 : t1_average1, name2 : t2_average1}) 
  avg2 = pd.DataFrame({name1 : t1_average2, name2 : t2_average2}) 
  ax1 = avg1.plot(kind='bar', legend=True, ax=axes[0]) 
  ax2 = avg2.plot(kind='bar', legend=False, ax=axes[1])
  axes[0].set_title('Computation Time')
  axes[1].set_title('Communication Time')
  ax1.set_ylabel('Time (seconds)')
  #ax2.set_xlabel('Time (seconds)')

  plt.tight_layout()
  fig.savefig(path)


def compare_max_by_section(times1, times2, path, name1='old', name2='new'):
  """Creates stack barplot of max computation and communication times 
  across all ranks.

  Parameters:
  -----------
  * times1 - Dict output from read_timer()
  * times2 - Dict output from read_times()
  * path - Path/name of output image
  """
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,10))

  # Calc average times first 
  sections1 = times1['comp_times'][COL_PLOT]
  sections2 = times1['comm_times'][COL_PLOT]
  t1_average1 = sections1.mean(0)
  t1_average2 = sections2.mean(0)

  sections1 = times2['comp_times'][COL_PLOT]
  sections2 = times2['comm_times'][COL_PLOT]
  t2_average1 = sections1.max(0)
  t2_average2 = sections2.max(0)

  # Plot
  avg1 = pd.DataFrame({name1 : t1_average1, name2 : t2_average1}) 
  avg2 = pd.DataFrame({name1 : t1_average2, name2 : t2_average2}) 
  ax1 = avg1.plot(kind='bar', legend=True, ax=axes[0]) 
  ax2 = avg2.plot(kind='bar', legend=False, ax=axes[1])
  axes[0].set_title('Computation Time')
  axes[1].set_title('Communication Time')
  ax1.set_ylabel('Time (seconds)')
  #ax2.set_xlabel('Time (seconds)')

  plt.tight_layout()
  fig.savefig(path)
