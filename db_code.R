#' Get the table information for a postgres database
#' http://tagteam.harvard.edu/hub_feeds/1981/feed_items/335812 source
#' @param config the configuration list
#' @return the table names, columns, and column types of all columns in the database
getTableInformation <- function(config = config.gp) {
  tables <- fetchQuery(
    "SELECT table_name, column_name, data_type 
    FROM information_schema.columns 
    WHERE table_name NOT LIKE '%prt%' 
    AND table_name NOT LIKE '%ext%' 
    AND table_name NOT LIKE '%tmp%' 
    ORDER BY 1, 2",
    config
  )
}

#' Replacement of the normal update function, you don't need to call this.
update <- function(object, ...) {
  args <- list(...)
  for (nm in names(args)) {
    object[[nm]] <- args[[nm]]
  }
  if (is.null(object$select)) {
    if (is.ident(object$from)) {
      var_names <- object$select
    }
    else {
      var_names <- qry_fields(object$src$con, object$from)
    }
    vars <- lapply(var_names, as.name)
    object$select <- vars
  }
  object$query <- dplyr:::build_query(object)
  object
}
#' Function to reflect a database, generalizable to others beyond postgres 
#' by simply changing getTableInformation appropriately
reflectDatabase <- function(config, envir.name = "tables",
                            subclass = "postgres") { # Note changed to sqlite subclass = "sqlite") { #
  if (!(envir.name %in% search())) {
    envir <- new.env(parent = .GlobalEnv)
  } else {
    envir <- as.environment(envir.name)
    detach(envir.name, character.only = TRUE)
  }
  src <- do.call(src_postgres, config)
  tables <- getTableInformation(config)
  tables <- split(tables, tables$table_name)
  lapply(tables, function(i) {
    nm <- ident(i$table_name[1])
    vars <- lapply(i$column_name, as.name)
    tbl <- dplyr::make_tbl(c(subclass, "sql"), src = src, from = nm,
                           select = vars, summarise = FALSE, mutate = FALSE,
                           where = NULL, group_by = NULL, order_by = NULL)
    tbl <- update(tbl)
    assign(
      nm,
      tbl,
      envir = envir
    )
  })
  attach(envir, name = envir.name)
}

searchTables <- function(str, env = "tables") {
  all.tbls <- ls(env)
  all.tbls[grep(str, all.tbls)]
}

reflectDatabase()
