#%%
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sns.set()

here = Path(".")
data = here / "data"
data.mkdir(exist_ok=True, parents=True)
#%% Part 1

# This is a dataframe (like tibble).
spotify_songs = pd.read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-01-21/spotify_songs.csv")
spotify_songs.to_csv(data / 'spotify_songs.csv', index=False)

spotify_songs
#%%
## How many songs are in each genre?
spotify_songs.groupby('playlist_genre').size()

#%%
## What is average value of energy and acousticness in the latin genre in this dataset?
spotify_songs[spotify_songs['playlist_genre'] == 'latin'][['energy', 'acousticness']].mean()

## Calculate the average duration of song (in minutes) across all subgenres. Which subgenre has the longest song on average?
spotify_songs.groupby('playlist_subgenre')['duration_ms'].mean().sort_values(ascending=False) / 60000

# Make two boxplots side-by-side of the danceability of songs stratifying by whether a song has a fast or slow tempo.
# Define fast tempo as any song that has a tempo above its median value.
median = spotify_songs['tempo'].median()
spotify_songs['tempo_type'] = np.where(spotify_songs['tempo'] > median, 'fast', 'slow')

sns.boxplot(spotify_songs, x='tempo_type', y='danceability')

#%%
sns.violinplot(spotify_songs, x='tempo_type', y='danceability')

#%% [markdown]
# ## Part 2
# Use this [tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html).
# The dataset can be procured as followed:
#%%
adata = sc.datasets.pbmc3k()
adata

# %%
