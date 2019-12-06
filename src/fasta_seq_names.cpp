/*
 *
 *	Author: Schuyler D. Smith
 *
 */

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <filesystem>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>


//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::List fasta_seq_names(
	std::string fasta_file_path
){
	Rcpp::StringVector seq_names;
	std::string 	seq_ID,
					not_seq_ID,
					line;

	Rcpp::Function directory("dir");
	Rcpp::Function file_path("file.path");
	Rcpp::Function basename("basename");
	Rcpp::Function file_path_sans_ext("file_path_sans_ext"); // must library(tools) in R code

	Rcpp::StringVector fasta_files = directory(fasta_file_path);
	Rcpp::StringVector namevec;
	if(fasta_files.length() > 0)
	{
		fasta_files = file_path(fasta_file_path, directory(fasta_file_path));
		namevec = file_path_sans_ext(directory(fasta_file_path));

	} else {
		fasta_files = file_path(fasta_file_path);
		namevec = file_path_sans_ext(basename(fasta_file_path));
	}

	size_t n_files = fasta_files.length();
	Rcpp::List fasta_seqs(n_files);
	for(size_t i=0; i < n_files; ++i)
	{
		std::ifstream 	fasta_file(fasta_file_path);
	    while(std::getline(fasta_file, line))
	    {
	    	std::stringstream   line_(line);
	    	std::getline(line_, seq_ID, '\t');
	    	line_ >> not_seq_ID; 
	    	if(seq_ID[0] == '>')
	    	{
	    		seq_names.push_back(seq_ID.erase(0, 1));
	    	}
	    }
    	fasta_seqs[i] = seq_names;
	}
	fasta_seqs.attr("names") = namevec;
return(fasta_seqs);
}