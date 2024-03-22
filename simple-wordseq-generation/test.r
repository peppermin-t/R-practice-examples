# split_func <- function(words, punc) {
#   punc_idx <- grep(punc, words, fixed = TRUE)
#   new_words <- rep("", length(words) + length(punc_idx))
#   words[punc_idx] <- gsub(punc, "", words[punc_idx], fixed = TRUE)

#   new_punc_idx <- punc_idx + 1:length(punc_idx)
#   new_words[new_punc_idx] <- punc
#   new_words[-new_punc_idx] <- words
#   new_words
# }

# a <- c("Hello,", "sdhjdf", "sdfsdfss,", "sdfsdfss,", "sdfs", "Hello,")
# # test for split_punct
# # print(split_func(a, ","))

# # test for unique
# print(unique(a))

# x <- cbind(c(1,2,1), c(2,2,2), c(3,4,4))
# x[c(FALSE, FALSE, FALSE), ]

# a <- c("apple,", "sdfsd,", "asdfasd")
# ipx <- grep(",", a)
# a[ipx]

# tabulate(c(3,4,4,6))

cat(c(1, 2))