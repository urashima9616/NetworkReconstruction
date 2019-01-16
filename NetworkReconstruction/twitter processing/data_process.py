import networkx as nx
import collections

#1. Parse all the tweets
#format:
#user_id 1
#follower 8
#following 9 
#is_retweet 20
#retweet_user_id 21

#user information
user_db = collections.defaultdict(list)

#user id to node idx and vice versa
id2idx = {}
idx2id = {}

#edge information
links = collections.defaultdict(list)

#degree information
retweeted_degree = collections.defaultdict(int)
max_degree = 0
max_candidate = None
def parse(tweet):
    i = 0 
    get_follower = 0
    follower = 0
    followee = 0
    retweet = 0
    retweet_id = ""
    while i < len(tweet):
        #check if it is first date
        if get_follower == 0 and len(tweet[i]) >= 12 and (tweet[i][1:4] == "200" or tweet[i][1:4] == "201") and "-" in tweet[i] :
            follower = int(tweet[i-2][1:len(tweet[i-2])-1])
            followee = int(tweet[i-1][1:len(tweet[i-1])-1])
            get_follower = 1
        if tweet[i] == '"false"' or tweet[i] == '"true"':
            retweet = 1 if tweet[i] == '"true"' else 0
            retweet_id = tweet[i+1][1:len(tweet[i+1])-1]
            break
        i += 1
    return follower, followee, retweet, retweet_id

idx = 0
count = 0
link_count = 0

with open("ira_tweets_csv_hashed.csv") as f:
    header = 0
    for line in f:
        if header == 0:
            header = 1
            continue
        else:
            tweet = line.split(",")
            user_id = tweet[1][1:len(tweet[0])-1]
            follower, followee, retweet, retweet_id = parse(tweet)
            #follower = int(users[8][1:len(users[8])-1])
            #followee = int(users[9][1:len(users[9])-1])
            #retweet =  users[20][1:len(users[20])-1]
            #retweet_id = users[21][1:len(users[21])-1]

            if user_id not in user_db:
                user_db[user_id].append(follower)
                user_db[user_id].append(followee)

            if user_id not in id2idx:
                id2idx[user_id] = idx
                idx2id[idx] = user_id
                idx += 1
            if retweet == 1:
                if retweet_id not in id2idx:
                    id2idx[retweet_id] = idx
                    idx2id[idx] = retweet_id
                    idx += 1
                links[retweet_id] = user_id
                retweeted_degree[retweet_id] += 1
                if retweeted_degree[retweet_id] > max_degree:
                    max_degree = retweeted_degree[retweet_id]
                    max_candidate = retweet_id
                link_count += 1
        if count%100000 == 0:
            print("Processed ", count, "tweets so far")
        count += 1

suspect_congr = {}
with open("exhibit_b.txt") as f:
    for each in f:
        lines = each.rstrip().split(" ")
        if lines:
            user_id = lines[0]
            if user_id in retweeted_degree:
                suspect_congr[user_id] = retweeted_degree[user_id]


print("Total ", len(user_db), " users recorded")
print("Total ", link_count, " links recorded")
print("Max degree is ", max_degree, " with user_id ", max_candidate)
print("suspect confirmed:  ", len(suspect_congr))

retwt_degree = sorted([(retweeted_degree[uid], uid) for uid in retweeted_degree.keys()])
with open("degree_distribution.csv", "w") as f:
    f.write("node_id,uid,degree\n")
    for each in retwt_degree:
        f.write(str(id2idx[each[1]])+","+each[1] + "," + str(each[0])+"\n")

with open("degree_distribution_suspect.csv", "w") as f:
    f.write("node_id, uid,degree\n")
    for each in retwt_degree:
        if each[1] in suspect_congr:
            f.write(str(id2idx[each[1]])+","+each[1] + "," + str(each[0])+"\n")







#2. Parse all the trolls

#3. Parse the tweets to construct the network
#user_id 1
#follower 7
#following 8 
#is_retweet 18
#retweet_user_id 19
