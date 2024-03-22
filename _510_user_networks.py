#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: 510_user_networks.py
---------------------------
This script reads hydrated tweet data from CSV files, creates a user network based
on mentions, and writes the network data to new CSV files. The network is
represented as an adjacency list, where each row represents an edge in the network.
The columns are the source user (i), the target user (j), and the weight of the
edge (w), which is the number of times i mentioned j in their tweets.

Author: Nurzhan Sapargali <nurzh.sapargali@gmail.com>
Date created: 2024-02-24
"""
from collections import Counter

import pandas as pd

HYDRATED_TWEETS = "./_900_output/data/hydrated/hydrated_tweets_{}.csv"
FILE_NUMBERS = list(range(250, 257)) + list(range(980, 987))
USER_NETS = "./_900_output/data/user_nets/user_net_{}.csv"


def main():
    # Loop over each file number
    for file_no in FILE_NUMBERS:
        print(file_no)

        # Initialize a Counter to store the edges of the network
        edges = Counter()

        # Read the hydrated tweet data, drop rows where the mentions is NaN
        df = pd.read_csv(HYDRATED_TWEETS.format(file_no), sep = "\t") 
        df = df.dropna(subset = ["mentions"])

        for i, row in df.iterrows():
            # Get the author of the tweet
            author = str(row["author_id"])
            # Get the list of users mentioned in the tweet from string of array
            mentions = row["mentions"][1:-1].replace("'", '').split(', ')
            # Create a list of edges from the author to each mentioned user
            neighbours = [(author, m) for m in mentions]
            # Update the Counter with the new edges
            edges.update(neighbours)
            
        # Create a new dataframe from the Counter
        adj_df = pd.DataFrame(
            [(i[0][0], i[0][1], i[1]) for i in edges.items()],
            columns = ["i", "j", "w"]
        )
        # Write the dataframe to a CSV file
        adj_df.to_csv(USER_NETS.format(file_no), index = False)


if __name__ == "__main__":
    main()
