#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: 500_hydrate_tweets.py
---------------------------
This script hydrates tweet IDs, i.e. it fetches the full data for each tweet
ID. It reads tweet IDs from CSV files named "corona_tweets_<c>.csv", where <c> is
a number. The data is fetched in batches of 100 IDs at a time, as this is the
maximum allowed by the Twitter API. If the rate limit is reached, it waits for 15
minutes (the standard rate limit window) before trying again. The hydrated tweet
data is saved to a new CSV file named "hydrated_tweets_<c>.csv

Author: Nurzhan Sapargali <nurzh.sapargali@gmail.com>
Date created: 2024-02-24
"""

import tweepy
import pandas as pd
from time import sleep

# Initialize the Tweepy client with your Twitter API credentials
client = tweepy.Client(consumer_key="",
                       consumer_secret="",
                       access_token="",
                       access_token_secret="")

for c in [250, 251, 252, 253, 254, 255, 256, 980, 981, 982, 983, 984, 985, 986]:
    dfs = []  # List to store dataframes

    # Read the raw tweet IDs
    raw = pd.read_csv(
        "./_200_input/tweets/corona_tweets_{}.csv".format(c),
        header=None
    )

    # Create indices for slicing the dataframe into batches of 100
    indices = list(range(0, raw.shape[0], 100))  
    for i in indices:
        # Slice the dataframe
        if i != indices[-1]:
            cut = raw.iloc[i:(i + 100), 0].to_list()
        else:
            cut = raw.iloc[i:, 0].to_list()

        timeout = True
        while timeout:
            try:
                # Fetch the tweet data from the Twitter API
                response = client.get_tweets(
                    cut,
                    tweet_fields = [
                        "author_id",
                        "in_reply_to_user_id",
                        "referenced_tweets"
                    ],
                    expansions = ["entities.mentions.username"],
                    user_auth = True
                )
                timeout = False

            except tweepy.TooManyRequests:
                # If we hit the rate limit, wait for 15 minutes
                sleep(15 * 64)
            
            except Exception as e:
                print(e)
                sleep(30)
        
        out = []
        # Process the response and append to the list
        for tweet in response.data:
            row = {}

            row["tweet_id"] = str(tweet.id)
            row["author_id"] = str(tweet.author_id)

            if tweet.in_reply_to_user_id:
                row["in_reply_to_user_id"] = str(tweet.in_reply_to_user_id)
            
            if tweet.referenced_tweets:
                for ref in tweet.referenced_tweets:
                    row[ref.type] = str(ref.id)

            if tweet.entities:
                ments = [m.get("id") for m in tweet.entities.get("mentions")]
                row["mentions"] = ments
            
            out.append(row)
        dfs.append(pd.DataFrame(out))
        print(len(dfs))
    
    # Concatenate all dataframes and save to a CSV file
    overall = pd.concat(dfs)
    overall.to_csv(
        "./_900_output/data/hydrated/hydrated_tweets_{}.csv".format(c),
        sep = "\t",
        index = False
    )
