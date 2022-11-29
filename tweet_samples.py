#!/usr/bin/env python3
# Authors:   Nurzhan Sapargali <nurzh.sapargali@gmail.com>

import tweepy
import pandas as pd

from time import sleep

client = tweepy.Client("")
dfs = []
raw = pd.read_csv("./_200_input/tweets/corona_tweets_982.csv", header=None)
indices = list(range(0, raw.shape[0], 100))
for i in indices:
    if i != indices[-1]:
        cut = raw.iloc[i:(i + 100), 0].to_list()
    else:
        cut = raw.iloc[i:, 0].to_list()
    timeout = True
    while timeout:
        try:
            response = client.get_tweets(cut,
                                         tweet_fields=["author_id", "in_reply_to_user_id", "referenced_tweets"],
                                         expansions=["entities.mentions.username"])
            timeout = False
        except Exception as e:
            print(e)
            if e == tweepy.TooManyRequests:
                sleep(15 * 64)
            else:
                sleep(30)
    out = []
    for tweet in response.data:
        row = {}
        row["tweet_id"] = str(tweet.id)
        #row["text"] = str(tweet.text)
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
overall = pd.concat(dfs)
overall.to_csv("./_900_output/data/hydrated/hydrated_tweets_982.csv",
               sep="\t",
               index=False)
    

# NETWORK CONSTRUCTION
# tot = pd.concat([pd.read_csv(i, sep="\t", lineterminator="\n") for i in glob("./_200_input/tweets/tweet_batch*.csv")])
# test = tot[~tot["mentions"].isna()]
# test["mentions"] = test["mentions"].apply(eval)
# decoup = pd.DataFrame(test["mentions"].values.tolist()).add_prefix("mention_")
# decoup.reset_index(drop=True, inplace=True)
# test.reset_index(drop=True, inplace=True)
# decoup["author_id"] = test["author_id"].astype(str)
# adj_lists = []
# for c in decoup.columns:
#     if "mention" in c:
#         cut = decoup[[c, "author_id"]]
#         cut = cut.dropna()
#         cut = cut[cut[c].isin(cut["author_id"])]
#         cut.columns = ["mention", "author"]
#         adj_lists.append(cut[["author", "mention"]])
