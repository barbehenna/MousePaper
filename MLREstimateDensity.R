simres = fread("~/Documents/Projects/MousePaper/data/trainingdata.csv")
simres = na.omit(simres)

num_train_examples = 70000
num_test_examples = 30000

train = simres[sample(x = nrow(simres), size = num_train_examples)]
test = simres[!(uuid %in% train$uuid)][sample(x=nrow(simres)-num_train_examples, size = num_test_examples)]

mod = lm(Density ~ square1 + square2 + square3 + square4 + square5 + square6 + square7 + square8, train)
summary(mod) # lots of colineraity expected

train_loss = mean((train$Density - fitted(mod))^2)
train_loss

test_loss = mean((test$Density - predict(mod, newdata = test))^2)
test_loss

