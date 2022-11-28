#!/usr/bin/env python3
# Authors:   Nurzhan Sapargali <nurzh.sapargali@gmail.com>

import tweepy
import pandas as pd

from time import sleep

client = tweepy.Client("")
#dfs = []
raw = pd.read_csv("./_200_input/tweets/corona_tweets_983.csv", header=None)
indices = list(range(450000, raw.shape[0], 100))
for i in indices:
    if i != indices[-1]:
        cut = raw.iloc[i:(i + 100), 0].to_list()
    else:
        cut = raw.iloc[i:, 0].to_list()
    timeout = True
    #if len(dfs) % 300 == 0:
        #sleep(15 * 60)
    while timeout:
        try:
            response = client.get_tweets(cut,
                                         tweet_fields=["author_id", "in_reply_to_user_id", "referenced_tweets"],
                                         expansions=["entities.mentions.username"])
            timeout = False
        #except tweepy.TooManyRequests:
        except Exception: 
            sleep(15 * 60)
    out = []
    for tweet in response.data:
        row = {}
        row["tweet_id"] = tweet.id
        row["text"] = tweet.text
        row["author_id"] = tweet.author_id
        row["in_reply_to_user_id"] = tweet.in_reply_to_user_id
        if tweet.referenced_tweets:
            for ref in tweet.referenced_tweets:
                row[ref.type] = ref.id
        if tweet.entities:
            ments = [m.get("id") for m in tweet.entities.get("mentions")]
            row["mentions"] = ments
        out.append(row)
    dfs.append(pd.DataFrame(out))
    print(len(dfs))