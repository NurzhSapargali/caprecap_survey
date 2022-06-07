#!/usr/bin/env python3
# Authors:   Nurzhan Sapargali <nurzh.sapargali@gmail.com>
"""Writes a 1% random sample of tweets to a file"""

class TweetSampler(tweepy.StreamingClient):
    def __init__(self, token):
        super(TweetSampler, self).__init__(token)
        self.data = []
    def on_data(self, tweets):
        try:
            raw = json.loads(tweets)
            self.data.append(raw)
        except BaseException as e:
            print(e)
            
stream = TweetSampler(key)
stream.sample(expansions="author_id")