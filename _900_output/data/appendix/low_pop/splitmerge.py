import pandas as pd
data = pd.read_csv("estimates_2.0_low_pop.csv")
nb = pd.read_csv("estimates_2.0_low_pop_mplenb.csv")
data = data[data["type"] != "MPLE-NB"]
pls = pd.concat([nb, data])
pls.head()
pls.to_csv("estimates_2.0_low_pop.csv", index=False)
