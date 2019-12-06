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
#include <RcppArmadillo.h>
#include <RcppParallel.h>

std::vector<std::string> extract_keys(std::unordered_map<std::string, double> const& input_map) 
{
	std::vector<std::string> values;
	for (auto const& element : input_map) 
	{
  		values.push_back(element.first);
  	}
  return values;
}

std::vector<double> extract_values(std::unordered_map<std::string, double> const& input_map) 
{
	std::vector<double> values;
	for (auto const& element : input_map) 
	{
		values.push_back(element.second);
	}
	return values;
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::DataFrame BLAST_primer_match(
	std::string BLAST_file_path,
	std::string primer_file_path,
	int min_length = 100,
	int min_id = 98
){

	std::unordered_map<std::string, double> primers;

	std::string 	query,
					prior_query,
					db_match,
					e_value,
					line,
					line_,
					primer_ID,
					not_primer_ID;;

	int 			length,
					mismatch,
					gap,
					qstart,
					qend,
					sstart,
					send,
					count;

	double 			perc_id,
					bitscore,
					prior_perc_id,
					prior_bitscore;

	std::ifstream 	primer_file(primer_file_path);
    while(std::getline(primer_file, line))
    {
    	std::stringstream   line_(line);
    	std::getline(line_, primer_ID, '\t');
    	line_ >> not_primer_ID; 
    	if(primer_ID[0] == '>')
    	{
    		primers[primer_ID.erase(0, 1)] = 0;
    	}
    }

	Rcpp::Function directory("dir");
	Rcpp::Function file_path("file.path");
	Rcpp::Function basename("basename");
	Rcpp::Function file_path_sans_ext("file_path_sans_ext"); // must library(tools) in R code

	Rcpp::StringVector BLAST_files = directory(BLAST_file_path);
	Rcpp::StringVector namevec;
	if(BLAST_files.length() > 0)
	{
		BLAST_files = file_path(BLAST_file_path, directory(BLAST_file_path));
		namevec = file_path_sans_ext(directory(BLAST_file_path));
	} else {
		BLAST_files = file_path(BLAST_file_path);
		namevec = file_path_sans_ext(basename(BLAST_file_path));
	}
	namevec.push_front("Primer");
	size_t n_files = BLAST_files.length();

	Rcpp::List primer_counts(1 + n_files);
	primer_counts[0] = extract_keys(primers);

	for(size_t i=0; i < n_files; ++i)
	{
		std::unordered_map<std::string, double> found_primers = primers;
		std::ifstream 	BLAST_file(BLAST_files[i]);
	    while(std::getline(BLAST_file, line)) 
	    {
	    	std::stringstream   line_(line);
	    	std::getline(line_, query, '\t');
	    	line_ >> db_match >> perc_id >> length >> mismatch >> gap >> 
	    		qstart >> qend >> sstart >> send >> e_value >> bitscore;
	    	if(query != prior_query & 
	    		length >= min_length &
	    		perc_id >= min_id)
	    	{
				if(found_primers.count(db_match))
				{
					found_primers[db_match]++;
				} 
	    	}
	    	prior_query = query;
	    }
    	primer_counts[i+1] = extract_values(found_primers);
	}
	primer_counts.attr("names") = namevec;

	Rcpp::DataFrame return_dataframe(primer_counts);
return(return_dataframe);
}

