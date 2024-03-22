# Authors: Yinjia Chen (s2520995@ed.ac.uk),
#          Yiwen Xing (s2530703@ed.ac.uk),
#          Chaoyue Yan (s2572326@ed.ac.uk)
# Contributions:
#   Joint completion: We worked on the first, second, and third problems
#       together, and worked on the comment section and the final code
#       simplification. We checked the validity of the code and made the
#       appropriate changes.
#   Yinjia Chen: Completed questions 8, 9 and 10 independently.
#   Yiwen Xing: Completed questions 6 and 7 independently.
#   Chaoyue Yan: Completed questions 4 and 5 independently.

# Overview: This programme generates a 50-word article based on the
#           given Ulysses text, under trigram/bigram/pure frequency
#           rules of word generation.

# setwd("txtfiles") ## comment out of submitted
a <- scan("4300-0.txt", what = "character", skip = 73, nlines = 32858 - 73)  ## read the text into a
a <- gsub("_(", "", a, fixed = TRUE)  ## remove "_("

split_punc <- function(words, punc) {
  # function to separate a specific punctuation from words in the word vector
  #     and put it behind
  # params: words: word vector to drag punctuations out
  #         punc:  the target punctuation to be pulled out
  punc_idx <- grep(punc, words, fixed = TRUE)  # find the index of words in "words" containing "punc"
  new_words <- rep("", length(words) + length(punc_idx))  # create a all-blank-strings vector, with size of the original a plus count of targeted punctuations
  words[punc_idx] <- gsub(punc, "", words[punc_idx], fixed = TRUE)  # remove "punc" in their original appearing location.
  # Note that fied=TRUE ensures punc to not be recognized as regular expressions

  new_punc_idx <- punc_idx + 1:length(punc_idx)  # specify the new single punc position in the new vector
  new_words[new_punc_idx] <- punc  # insert the punc in the right position
  new_words[-new_punc_idx] <- words  # insert the words in the corresponding else where
  return(new_words)  # return the new word vector
}

punc_to_split <- c(",", ".", ";", "!", ":", "?")

for (punc in punc_to_split) {
  a <- split_punc(a, punc)  # pull punc out with defined function split_punc
}
# obtain the final word vector, a, with punctuations pulled out

a_lower <- tolower(a)  # a with all words in lower case
dict <- unique(a_lower)  # dictionary of the article, i.e. each unique word only occurs once here
a_idx <- match(a_lower, dict)  # a_lower with each word represented as its word order in dict
freq <- tabulate(a_idx)  # the occuring times of each word in the a_lower, in the word order of dict

thr <- 0
repeat {  # searching process of the threshold of sequence of "common words"
  # seach starts from 0
  common_word_counts <- length(which(freq > thr))  # obtain the the count of common words with current threshold
  if (common_word_counts < 1000) {
    # if the common word count first exceed 1000, judge:
    if (1000 - common_word_counts > length(which(freq > thr - 1)) - 1000) {  
      # thr & thr - 1 are two thresholds getting common word count around 1000 most closely, compare which one is closer
      thr <- thr - 1  # if thr - 1 is closer, set thr - 1
    }
    break
  }
  thr <- thr + 1  # iterate thr
}
b <- dict[which(freq > thr)]  # common word dictionary with threshold thr (different from "dict")

common_words_idx <- match(a_lower, b)  # a_lower with each word represented as its word order in b (display NA for words not in b)

col1 <- common_words_idx  # first column in T
col2 <- c(common_words_idx[-1], NA)  # second column in T (shift up one step)
col3 <- c(common_words_idx[-c(1, 2)], c(NA, NA))  # third column in T (shift up teo steps)
T <- cbind(col1, col2, col3)  # combine three columns into one matrix (each row: three adjacent words in article, i.e. a word triplet in article)
# cbind returns a matrix, not df
T <- T[which(!is.na(rowSums(T))), ]  # delete rows containing NA

P <- cbind(col1, col2)  # combine two columns into one matrix (each row: two adjacent words in article, i.e. a word pair in article)
P <- P[which(!is.na(rowSums(P))), ]  # delete rows containing NA

generate <- function(n, T, P, common_words_idx) {
  # function to generate an n-words article based on the given
  #     texts under tri-gram, bigram, word frequency rules
  # params: n:                count of generating words
  #         T:                the trigram matrix of the article text
  #         P:                the bigram matrix of the article text
  #         common_words_idx: article in lower case with each word
  #                           represented as its word order in b
  last2 <- NA  # the second previous word
  last1 <- NA  # the previous word
  article <- rep(NA, n)  # n-legth vector for storing the newly generated article
  for (i in 1:n) {
    tri_gram_candidates <-
      matrix(T[which(T[, 1] == last2 & T[, 2] == last1), ], ncol = 3)  # subtract the triplet submatrix with the words in first two columns equal to the given previous two words
    bi_gram_candidates <-
      matrix(P[which(P[, 1] == last1), ], ncol = 2)  # subtract the pair submatrix with the words in first two columns equal to the given previous two words
    # matrix(): demanding bi & tri to be 2-d matricies even if having no rows
    if (nrow(tri_gram_candidates) != 0) {
      candidates <- tri_gram_candidates[, 3]  # if rows exist in triplet submatrix, subtract candidate words in the third column
    } else if (nrow(bi_gram_candidates) != 0) {
      candidates <- bi_gram_candidates[, 2]  # else if rows exist in pair submatrix, subtract candidate words in the second column
    } else {
      candidates <- common_words_idx[which(!is.na(common_words_idx))]  # else candidate words are all non-NA words
    }
    if (length(candidates) == 1) article[i] <- candidates
    else article[i] <- sample(x = candidates, size = 1, replace = TRUE)
    # sample a word in candidates with each of them in identical probability, thus the prob of picking a specific word is exactly the occuring times of that word in candidates
    # Note that replce=True means the chosen elsement can be chosen again next time
    # concatenate the selected word in the new article
    last2 <- last1
    last1 <- article[i]  # renew the previous words for looping
  }
  cat(paste(b[article], collapse = " "), "\n")  # return the new article in a single string by joining each ward with a space
}

generate(50, T, P, common_words_idx)  # generate a 50-word article and print it, under the standard criteria

generate_only_from_freq <- function(n, common_words_idx) {
  # function to generate an n-words article based on the given texts
  #     solely under word frequencies
  # params: n:                count of generating words
  #         common_words_idx: article in lower case with each word
  #                           represented as its word order in b
  article <- rep(NA, n)
  for (i in 1:n) {
    candidates <- common_words_idx[which(!is.na(common_words_idx))]
    if (length(candidates) == 1) article[i] <- candidates
    else article[i] <- sample(x = candidates, size = 1, replace = TRUE)
  }
  cat(paste(b[article], collapse = " "), "\n")
}

generate_only_from_freq(50, common_words_idx)  # generate a 50-word article and print it, only based on word frequencies

common_words_cap_freq <- tabulate(match(a_lower[a != a_lower], dict))  # the occuring times of each capitalized word in the article, in the word order of dict
common_words_freq <- tabulate(match(a_lower, dict))  # the occuring times of each word in the a_lower, in the word order of dict

common_words_cap_prob <- common_words_cap_freq / common_words_freq  # the capitalizing prob of each word in the article, in the word order of dict
cap_word_set <- dict[which(common_words_cap_prob > 0.5)]  # the words with capitalizing prob greater than 0.5

capitalize <- function(word) {
  # function to capitalize a word
  # params: word: word to capitalize
  paste(toupper(substr(word, 1, 1)), substr(word, 2, nchar(word)), sep = "")
}

generate_cap <- function(n, T, P, common_words_idx) {
  # function to generate an n-words article based on the given
  #     texts under tri-gram, bigram, word frequency rules,
  #     with a specific threshold of capitalization probability
  # params: n:                count of generating words
  #         T:                the trigram matrix of the article text
  #         P:                the bigram matrix of the article text
  #         common_words_idx: article in lower case with each word
  #                           represented as its word order in b
  last2 <- NA
  last1 <- NA
  article <- rep(NA, n)
  for (i in 1:n) {
    tri_gram_candidates <-
      matrix(T[which(T[, 1] == last2 & T[, 2] == last1), ], ncol = 3)
    bi_gram_candidates <-
      matrix(P[which(P[, 1] == last1), ], ncol = 2)
    if (nrow(tri_gram_candidates) != 0) {
      candidates <- c(tri_gram_candidates[, 3])
    } else if (nrow(bi_gram_candidates) != 0) {
      candidates <- c(bi_gram_candidates[, 2])
    } else {
      candidates <- c(common_words_idx[which(!is.na(common_words_idx))])
    }
    if (length(candidates) == 1) article[i] <- candidates
    else article[i] <- sample(x = candidates, size = 1, replace = TRUE)
    last2 <- last1
    last1 <- article[i]
  }
  article <- b[article]
  cap_id <- which(article %in% cap_word_set)
  article[cap_id] <- capitalize(article[cap_id])
  cat(paste(article, collapse = " "), "\n")
}

generate_cap(50, T, P, common_words_idx)  # generate a 50-word article and print it, under the standard criteria with capitalization

# Instructor's comments:
# Would be better to put overview comments describing what the code is about near the start of the code file
# Don’t put question numbers in comments. Code should be understandable without the practical sheet.
# Functions should start with a comment describing inputs, purpose and outputs.
# better to put split_punc call in a loop
# It is really bad practice for functions to use data not supplied in the argument list. This is a really common way to write non-portable code that does not work as intended!
# Your 'generate' code is wrong because candidates can end up being a single integer so that sample(candidates,1) produces a number between 1 and candidates, which is not what is wanted.
# Your 'generate_cap' code appears to be further wrong as you are mixing up words and word indices. e.g. (ma_t[which(ma_t[, 1] == last2 & ma_t[, 2] == last1), ] - here ma_t contains integers and 'last1' and 'last2' character strings, so it can not be doing what is intended.
# Generally it would have been much better to generate a sequence of word indices, say ii, and then produce the equivalent words with something like b[ii].