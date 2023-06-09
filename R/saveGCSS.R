#' Save CSS File Embedded in Shiny App
#'
#' @inheritParams saveFiles
#' @returns No return value. It is used to save necessary files temporarily
#' to run Shiny App.

saveGCSS <- function(filePath){

  globstr <- "html, body {
	position: relative;
	width: 100%;
	height: 100%;
}

body {
	color: #333;
	margin: 0;
	box-sizing: border-box;
	font-family: 'Montserrat', 'Helvetica Neue', sans-serif;
}

a {
	color: rgb(0,100,200);
	text-decoration: none;
}

a:hover {
	text-decoration: underline;
}

a:visited {
	color: rgb(0,80,160);
}

label {
	display: block;
}

input, button, select, textarea {
	font-family: inherit;
	font-size: inherit;
	padding: 0.4em;
	margin: 0 0 0.5em 0;
	box-sizing: border-box;
	border: 1px solid #ccc;
	border-radius: 2px;
}

input:disabled {
	color: #ccc;
}

input[type='range'] {
	height: 0;
}

button {
	color: #333;
	background-color: #f4f4f4;
	outline: none;
}

button:active {
	background-color: #ddd;
}

button:focus {
	border-color: #666;
}
"
globsplt <- strsplit(globstr, "\n")
globvec <- unlist(globsplt)

fileConn <- file(paste(filePath, "/shinyapp/www/global.css", sep=""))
writeLines(globvec, fileConn)
close(fileConn)

}
