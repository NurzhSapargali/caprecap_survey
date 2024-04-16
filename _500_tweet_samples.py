#!/usr/bin/env python3
# Authors:   Nurzhan Sapargali <nurzh.sapargali@gmail.com>
import config

import tweepy
import pandas as pd

from time import sleep

client = tweepy.Client(
    consumer_key = config.CONSUMER_KEY,
    consumer_secret = config.CONSUMER_SECRET,
    access_token = config.ACCESS_TOKEN,
    access_token_secret = config.ACCESS_TOKEN_SECRET
)

for c in [251, 252, 253, 254, 255]:
    dfs = []
    raw = pd.read_csv("./_200_input/tweets/corona_tweets_{}.csv".format(c), header=None)
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
                                             expansions=["entities.mentions.username"],
                                             user_auth=True)
                timeout = False
            except tweepy.TooManyRequests:
                sleep(15 * 64)
            except Exception as e:
                print(e)
                sleep(30)
        out = []
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
    overall = pd.concat(dfs)
    overall.to_csv("./_900_output/data/hydrated/hydrated_tweets_{}.csv".format(c),
                   sep="\t",
                   index=False)
    

# NETWORK CONSTRUCTION
tot = [pd.read_csv(i, sep="\t", lineterminator="\n") for i in glob("./_900_output/data/hydrated/hydrated_tweets_9*.csv")]
clean = []
for n in tot:
  test = n[~n["mentions"].isna()]
  test["mentions"] = test["mentions"].apply(eval)
  decoup = pd.DataFrame(test["mentions"].values.tolist()).add_prefix("mention_")
  decoup.reset_index(drop=True, inplace=True)
  test.reset_index(drop=True, inplace=True)
  decoup["author_id"] = test["author_id"].astype(str)
  adj_lists = []
  for c in decoup.columns:
      if "mention" in c:
          cut = decoup[[c, "author_id"]]
          ser = cut.dropna()
          ser = ser[ser[c].isin(cut["author_id"])]
          ser.columns = ["mention", "author"]
          adj_lists.append(ser[["author", "mention"]])
  clean.append(pd.concat(adj_lists))
