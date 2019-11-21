# For processing we used mufold server:
# http://dslsrv2.eecs.missouri.edu/~zlht3/ss

# sequence at the time (~1m:30s processing time).

# get table of statuses and jobid if task is new than return empty table
get_status_all <- function(email){

  # Get page with results back using email identifier
  retrive_results_status <- function(email){
    httr::POST("http://dslsrv2.eecs.missouri.edu/~zlht3/ss/retrieve",
               body = list(jobid="",email=email),
               encode = c("form")) # , httr::verbose()

  }

  ## get results page with the table of finished sequence
  results <- retrive_results_status(email)

  ## extract table from results using xpath
  jobid_entries <- results %>% httr::content() %>%
    rvest::html_nodes(xpath='//html/body/div[2]/div[1]/div[2]/div/div/div[2]/table') %>%
    rvest::html_table()

  # add right header to the table
  jobid_entries <- jobid_entries[[1]]
  names(jobid_entries)<-c("num", "status", "target","time","len","sstype", "jobid")
  jobid_entries
}

# send sequence to server
send_to_server <- function(seq, email){
  lbody = list(email=email, target_name = "", 'tasks[]'=c(1), sequence=seq)
  httr::POST("http://dslsrv2.eecs.missouri.edu/~zlht3/ss/submit",
             body = lbody,
             encode = c("form"))

  # This delay not required, but we need be polite
  Sys.sleep(0.5)
}

# convert data.frame of sequences to vector of fasta lines by 10 (server limit on one upload)
df_to_fastaline <- function(df_sequences) {
  df_sequences %>%
    dplyr::rowwise() %>%
    dplyr::mutate(rec=paste0(">",target,"\n",fasta)) %>%
    # no more than 10 sequences at once
    dplyr::mutate(num = dplyr::row_number()%%round(1+dplyr::n()/10)) %>%
    dplyr::select(rec,num) %>%
    dplyr::group_by(num) %>%
    dplyr::summarise(fasta_line=stringr::str_c(c(rec), sep = "\n", collapse = "\n"))
}

#' Title Send sequences to server
#'
#' @param df_sequences data.frame with sequences (columns target and fasta required)
#' @param email at the moment any email-like identifier, it is not checked
#'
#' @return object which you can use as input to retrieve results from server
#' @export
#'
#' @examples
#' email = "test11@test.com"
#' testseq <- data.frame(target = c("fasta_seq1", "fasta_seq2"),
#'                       fasta  = c("RDPHPAPPRTSQEPHQNPHGVIPSESKPFVASKPKPHT",
#'                                  "PSLPRPPGCPAQQGGSAPIDPPPVHESPHPPLPATEPASRLSSE"),
#'                       stringsAsFactors=FALSE)
#'
#' # send sequences to server
#'
#' sequences <- sendToServer(testseq, email)
#'
#' # wait a couple of minutes, 1 sequence ~ 1-2 minutes of calculation
#'
#' # retrive results from server ( only finished sequences )
#'
#' structures <- getResultsFromServer(sequences)
#'
#' # or you can obtain all finished sequences which associated with this email
#'
#' all.structures <- getAllResultsFromServer(email)
sendToServer <- function (df_sequences, email) {

  # get current status and num
  accountStatus <- get_status_all(email)

  # next number for task
  nextNumberOfTask <- nrow(accountStatus)+1

  # vector of sequences
  preparedToSend <- df_to_fastaline(df_sequences)

  # send to server using map
  purrr::map2(preparedToSend$fasta_line, email,
              function(fasta_line, email) send_to_server(fasta_line, email))

  # return new object (list) with n, email and list of sequences
  result <- list(sequences = df_sequences)
  attr(result,'current_task') <- nextNumberOfTask
  attr(result,'email') <- email
  result
}

#' Get results from server MUFOLD server
#'
#' @param send_object object which was obtained by function 'sendToServer'
#'
#' @return data.frame with prediction of secondary structure
#' @export
#'
#' @examples
#' email = "test11@test.com"
#' testseq <- data.frame(target = c("fasta_seq1", "fasta_seq2"),
#'                       fasta  = c("RDPHPAPPRTSQEPHQNPHGVIPSESKPFVASKPKPHT",
#'                                  "PSLPRPPGCPAQQGGSAPIDPPPVHESPHPPLPATEPASRLSSE"),
#'                       stringsAsFactors=FALSE)
#'
#' # send sequences to server
#'
#' sequences <- sendToServer(testseq, email)
#'
#' # wait a couple of minutes, 1 sequence ~ 1-2 minutes of calculation
#'
#' # retrive results from server ( only finished sequences )
#'
#' structures <- getResultsFromServer(sequences)
#'
#' # or you can obtain all finished sequences which associated with this email
#'
#' all.structures <- getAllResultsFromServer(email)
getResultsFromServer <- function(send_object) {

  ### get Q3 and Q8 from remote server
  # http://dslsrv2.eecs.missouri.edu/~zlht3/ss/download_ss_results_only/ss_5d8936a4a0d86
  get_prediction <- function(jobid){
    url <- "http://dslsrv2.eecs.missouri.edu/~zlht3/ss/download_ss_results_only/"
    ss <- httr::GET(paste0(url, jobid))
    ss %>%
      httr::content(as="text", encoding = "UTF-8") %>%
      stringr::str_trim(.) %>%
      stringr::str_split(., "\n", simplify = T) %>%
      list(Q3 = .[,2], Q8=.[,4])
  }

  # safe get_prediction
  safe_get_prediction <- purrr::possibly(get_prediction,list(Q3=NA,Q8=NA), quiet = T)

  email <- attr(send_object,'email')
  current_task <- attr(send_object,'current_task')

  # get statuses
  tableOfResults <- get_status_all(email)

  # do we have still running jobs?
  runningTasks <- dplyr::filter(tableOfResults, num >= current_task, status != "Finished")
  if (nrow(runningTasks) != 0) {
    allTasks <- nrow(dplyr::filter(tableOfResults, num >= current_task))
    warning(paste0("\nWe still have running tasks: ",nrow(runningTasks),"/",allTasks))
  }

  # filter results
  finishedTasks <- dplyr::filter(tableOfResults, num >= current_task, status == "Finished")
  if (nrow(finishedTasks) > 0) {
    finishedTasks$target <- stringr::str_sub(finishedTasks$target, end=-3)

    # retrive structure prediction
    finishedStructures <- finishedTasks %>%
      dplyr::select(target, jobid) %>%
      purrr::pmap_df(., function(target, jobid) {
        Sys.sleep(0.3)
        print(jobid)
        lst <- safe_get_prediction(jobid)
        tibble::tibble(
          target = target,
          Q3 = lst$Q3,
          Q8 = lst$Q8
        )})

    dplyr::full_join(send_object$sequences, finishedStructures, by=c("target"))
  }
}

getAllResultsFromServer <- function(email) {
  fake_object <- list(sequences=data.frame(target=c(character()),
                                           fasta=c(character()),
                                           stringsAsFactors = F))
  attr(fake_object,"current_task") <- 1
  attr(fake_object, "email") <- email

  getResultsFromServer(fake_object)
}


# email = "testtest1@test.com"
# testseq <- data.frame(target = c("new_testqq",
#                                  "new_testqq2"),
#                       fasta = c("RDPHPAPPRTSQEPHQNPHGVIPSESKPFVASKPKPHT",
#                                 "PSLPRPPGCPAQQGGSAPIDPPPVHESPHPPLPATEPASRLSSE"),
#                       stringsAsFactors=FALSE)
#
# df_to_fastaline(testseq)
#
# strObj <- sendToServer(testseq, email)
# str(strObj)
#
# res <- getResultsFromServer(strObj)
# res
#
# res1 <- getAllResultsFromServer(email)
# res1


