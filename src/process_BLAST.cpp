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
Rcpp::DataFrame process_BLAST(
    std::string BLAST_file_path,
    int min_length = 100,
    int id_cutoff = 98
){
	std::string 	query,
					subject,
					e_value,
					line,
					line_,
					prior_query, 
					prior_subject,
					prior_e_value;

	int 			length,
					mismatch,
					gap,
					qstart,
					qend,
					sstart,
					send,
					count,
					prior_length,
					prior_mismatch,
					prior_gap,
					prior_qstart,
					prior_qend,
					prior_sstart,
					prior_send;

	double 			perc_id,
					bitscore,
					prior_perc_id,
					prior_bitscore;

	Rcpp::Function directory("dir");
	Rcpp::Function file_path("file.path");
	Rcpp::Function file_path_sans_ext("file_path_sans_ext"); // must library(tools) in R code

	Rcpp::CharacterVector BLAST_files = file_path(BLAST_file_path, directory(BLAST_file_path));

	size_t n_files = BLAST_files.length();
	std::vector<int> unique_reads(n_files);
	std::vector<int> multiple_alignments(n_files);

	std::vector<double> matches;
	std::vector<double> hundreds(n_files);
	std::vector<double> ninetynines(n_files);
	std::vector<double> ninetyeights(n_files);
	std::vector<double> ninetysevens(n_files);
	std::vector<double> ninetysixes(n_files);
	std::vector<double> ninetyfives(n_files);
	std::vector<double> ninetythrees(n_files);
	std::vector<double> ninetys(n_files);
	std::vector<double> eightyfives(n_files);
	std::vector<double> eightys(n_files);
	std::vector<double> seventys(n_files);
	std::vector<double> sixtys(n_files);

	for(size_t i=0; i < n_files; ++i)
	{
		int unique_read = 0;
		int multiple_alignment = 0;
		std::vector<double> read_perc_id;

		std::ifstream 	BLAST_file(BLAST_files[i]);
	    while(std::getline(BLAST_file, line)) 
	    {
	    	std::stringstream   line_(line);
	    	std::getline(line_, query, '\t');
	    	line_ >> subject >> perc_id >> length >> mismatch >> gap >> 
	    		qstart >> qend >> sstart >> send >> e_value >> bitscore;
	    	if(length >= min_length)
		    {
			    if(query != prior_query & perc_id >= 60)
		   		{
		    		++unique_read;
		    		read_perc_id.push_back(perc_id);
		    	}
			    if(query == prior_query & perc_id >= id_cutoff & length == prior_length)
		   		{
		    		++multiple_alignment;
		    	}
	    	}

	    	prior_query = query;
	    	prior_subject = subject;
	    	prior_perc_id = perc_id;
	    	prior_length = length;
	    	prior_mismatch = mismatch;
	    	prior_gap = gap;
	    	prior_qstart = qstart; 
	    	prior_qend = qend;
	    	prior_sstart = sstart;
	    	prior_send = send;
	    	prior_e_value = e_value;
	    	prior_bitscore = bitscore;
	    }
	    unique_reads[i] = unique_read;
	    multiple_alignments[i] = multiple_alignment;

		hundreds[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									 [](double i){return i == 100;}
								);
		ninetynines[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 99);
								});
		ninetyeights[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 98);
								});
		ninetysevens[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 97);
								});
		ninetysixes[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 96);
								});
		ninetyfives[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 95);
								});
		ninetythrees[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 93);
								});
		ninetys[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 90);
								});
		eightyfives[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 85);
								});
		eightys[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 80);
								});
		seventys[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 70);
								});
		sixtys[i] = std::count_if(read_perc_id.begin(), 
									read_perc_id.end(), 
									[&](double const &i) {
										return (i >= 60);
								});

	}
return Rcpp::DataFrame::create(
    Rcpp::Named("File") = file_path_sans_ext(directory(BLAST_file_path)),
    Rcpp::Named("Unique_Reads") = unique_reads,
    Rcpp::Named("Multiple_Alignments") = multiple_alignments,
    Rcpp::Named("id_100") = hundreds,
    Rcpp::Named("id_99") = ninetynines,
    Rcpp::Named("id_98") = ninetyeights,
    Rcpp::Named("id_97") = ninetysevens,
    Rcpp::Named("id_96") = ninetysixes,
    Rcpp::Named("id_95") = ninetyfives,
    Rcpp::Named("id_93") = ninetythrees,
    Rcpp::Named("id_90") = ninetys,
    Rcpp::Named("id_85") = eightyfives,
    Rcpp::Named("id_80") = eightys,
    Rcpp::Named("id_70") = seventys,
    Rcpp::Named("id_60") = sixtys
    );
}