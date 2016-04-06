#include <vector>
#include <string>

#include "alignment_utilities.h"



std::string reverse(const std::string& str) { 
  std::string result = "";
  for (int i = 0; i < str.length(); ++i) 
    result = str[i] + result; 
  return result; 
}

// find what index to partition the string at
int partition(const std::vector<int>& ScoreL, const std::vector<int>& ScoreR) {
  int index = 0;
  int max_score = ScoreL[0] + ScoreR[ScoreR.size()-1];
  for (int i = 1; i < ScoreL.size(); ++i) {
    if (ScoreL[i] + ScoreR[ScoreR.size() - i - 1] >= max_score) {
      max_score = ScoreL[i] + ScoreR[ScoreR.size() - i - 1];
      index = i;
    }
  }
  return index;
}

// Get the sequence from the provided fasta file
bool getSequence(const std::string& file_name, std::string& seq) {
  std::ifstream sequenceInput(file_name.c_str());
  if (!sequenceInput) {
    std::cout << "Unable to open sequence file " << file_name << std::endl;
    return false;
  }
  std::string seq_in;
  getline(sequenceInput, seq_in); // discard the line of description
  seq_in.clear();
  while (sequenceInput >> seq_in)
  	seq += seq_in;
  sequenceInput.close();
  return true;
}

std::vector<int> NWScore(const std::string& seq1, const std::string& seq2) {
  int WordScore[2][seq2.length()+1];
  int gap = -1;  //gap penalty

  for (int i=0; i < 2; ++i)
    WordScore[i][0] = i*gap;

  for (int j=0; j < seq2.length()+1; ++j)
    WordScore[0][j] = j*gap;

  for (int i=1; i < seq1.length()+1; ++i){
    WordScore[i%2][0] = i*gap;
    for (int j=1; j < seq2.length()+1; ++j){
      //score for comparing the letters with each other
      int match = WordScore[(i-1)%2][j-1] + Score(seq1[i-1], seq2[j-1]);
      //score for comparing the letter in seq2 with a gap in seq1
      int gap_seq1 = WordScore[i%2][j-1] + gap; 
      //score for comparing the letter in seq1 with a gap in seq2
      int gap_seq2 = WordScore[(i-1)%2][j] + gap;
      //maximum score of the 2 possible gaps
      int maxGapScore = std::max(gap_seq1,gap_seq2);
      WordScore[i%2][j] = std::max(match,maxGapScore);//total maximum
    }
  }

  std::vector<int> last_line;
  for (int i = 0; i < seq2.length()+1; ++i) 
    last_line.push_back(WordScore[seq1.length()%2][i]);

  return last_line;
}

std::vector<std::string> NWAlignment(const std::string& seq1, const std::string& seq2) {
  int WordScore[seq1.length()+1][seq2.length()+1];
  int gap = -1;  //gap penalty

  for (int i=0; i < seq1.length()+1; ++i)
    WordScore[i][0] = i*gap;

  for (int j=0; j < seq2.length()+1; ++j)
    WordScore[0][j] = j*gap;

  for (int i=1; i < seq1.length()+1; ++i){
    for (int j=1; j < seq2.length()+1; ++j){
      //score for comparing the letters with each other
      int match = WordScore[i-1][j-1] + Score(seq1[i-1], seq2[j-1]);
      //score for comparing the letter in seq2 with a gap in seq1
      int gap_seq1 = WordScore[i][j-1] + gap; 
      //score for comparing the letter in seq1 with a gap in seq2
      int gap_seq2 = WordScore[i-1][j] + gap;
      //maximum score
      WordScore[i][j] = std::max(match, std::max(gap_seq1,gap_seq2));//total maximum
    }
  }

  std::vector<std::string> alignment = {"","",""};
  int i = seq1.length();
  int j = seq2.length();
  while (i > 0 or j > 0) {

    if (i > 0 and j > 0 and WordScore[i][j] == WordScore[i-1][j-1] + Score(seq1[i-1], seq2[j-1])) {
      alignment[0] = seq1[i-1] + alignment[0];
      alignment[1] = (seq1[i-1] == seq2[j-1]) ? "|" + alignment[1] : " " + alignment[1];
      alignment[2] = seq2[j-1] + alignment[2];
      --i;
      --j;
    }

    else if (i > 0 and WordScore[i][j] == WordScore[i-1][j] + gap) {
      alignment[0] = seq1[i-1] + alignment[0];
      alignment[1] = " " + alignment[1];
      alignment[2] = "-" + alignment[2];
      --i;
    }

    else {
      alignment[0] = "-" + alignment[0];
      alignment[1] = " " + alignment[1];
      alignment[2] = seq2[j-1] + alignment[2];
      --j;
    }

  }
  return alignment;
}

std::vector<std::string> hirschbergAlignment(const std::string& seq1, const std::string& seq2) {
  std::vector<std::string> alignment = {"","",""};

  if (seq1.length() == 0) {
    for (int i = 0; i < seq2.length(); ++i) {
      alignment[0] += '-';
      alignment[1] += " ";
      alignment[2] += seq2[i];
    }
  }

  else if (seq2.length() == 0) {
    for (int i = 0; i < seq1.length(); ++i) {
      alignment[0] += seq1[i];
      alignment[1] += " ";
      alignment[2] += '-';
    }
  }

  else if (seq1.length() == 1 or seq2.length() == 1)
    alignment = NWAlignment(seq1,seq2);

  else {
    std::vector<int> ScoreL = NWScore(seq1.substr(0, seq1.length()/2), seq2);
    std::vector<int> ScoreR = NWScore(reverse(seq1.substr(seq1.length()/2)), reverse(seq2));

    int split_seq2_at = partition(ScoreL, ScoreR);

    std::vector<std::string> AlignL = hirschbergAlignment(seq1.substr(0, seq1.length()/2), seq2.substr(0,split_seq2_at));
    std::vector<std::string> AlignR = hirschbergAlignment(seq1.substr(seq1.length()/2), seq2.substr(split_seq2_at));
    alignment[0] = AlignL[0] + AlignR[0];
    alignment[1] = AlignL[1] + AlignR[1];
    alignment[2] = AlignL[2] + AlignR[2];
  }
  return alignment;
}