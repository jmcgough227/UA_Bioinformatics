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
  int gap = -2;  //gap penalty

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
  int gap = -2;  //gap penalty

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

/********************************************************************************/
void operator+=(alignmentStats &a, const alignmentStats &b) {
	a.insertion += b.insertion;
	a.deletion += b.deletion;
	a.synMutation += b.synMutation;
	a.nonSynMutation += b.nonSynMutation;
	a.match += b.match;
	a.mismatch += b.mismatch;
}


bool aminoCheck ( const std::string& str1, const std::string& str2 ) {
  std::map<std::string, std::string> RNA_codon_table = {
  // Second Base
  // U C A G
  // U
  {"UUU", "Phe"}, {"UCU", "Ser"}, {"UAU", "Tyr"}, {"UGU", "Cys"}, // UxU
  {"UUC", "Phe"}, {"UCC", "Ser"}, {"UAC", "Tyr"}, {"UGC", "Cys"}, // UxC
  {"UUA", "Leu"}, {"UCA", "Ser"}, {"UAA", "---"}, {"UGA", "---"}, // UxA
  {"UUG", "Leu"}, {"UCG", "Ser"}, {"UAG", "---"}, {"UGG", "Trp"}, // UxG
  // C
  {"CUU", "Leu"}, {"CCU", "Pro"}, {"CAU", "His"}, {"CGU", "Arg"}, // CxU
  {"CUC", "Leu"}, {"CCC", "Pro"}, {"CAC", "His"}, {"CGC", "Arg"}, // CxC
  {"CUA", "Leu"}, {"CCA", "Pro"}, {"CAA", "Gln"}, {"CGA", "Arg"}, // CxA
  {"CUG", "Leu"}, {"CCG", "Pro"}, {"CAG", "Gln"}, {"CGG", "Arg"}, // CxG
  // A
  {"AUU", "Ile"}, {"ACU", "Thr"}, {"AAU", "Asn"}, {"AGU", "Ser"}, // AxU
  {"AUC", "Ile"}, {"ACC", "Thr"}, {"AAC", "Asn"}, {"AGC", "Ser"}, // AxC
  {"AUA", "Ile"}, {"ACA", "Thr"}, {"AAA", "Lys"}, {"AGA", "Arg"}, // AxA
  {"AUG", "Met"}, {"ACG", "Thr"}, {"AAG", "Lys"}, {"AGG", "Arg"}, // AxG
  // G
  {"GUU", "Val"}, {"GCU", "Ala"}, {"GAU", "Asp"}, {"GGU", "Gly"}, // GxU
  {"GUC", "Val"}, {"GCC", "Ala"}, {"GAC", "Asp"}, {"GGC", "Gly"}, // GxC
  {"GUA", "Val"}, {"GCA", "Ala"}, {"GAA", "Glu"}, {"GGA", "Gly"}, // GxA
  {"GUG", "Val"}, {"GCG", "Ala"}, {"GAG", "Glu"}, {"GGG", "Gly"}  // GxG
  };
  return RNA_codon_table[str1] == RNA_codon_table[str2];
}


alignmentStats gatherStats( const std::vector<std::string>& alignment ) {
  
  alignmentStats data;
  
  std::string aminoStr0;
  std::string aminoStr1;
  std::string aminoStr2;
  
  bool gapFlag;
  
  for ( int i = 0; i < alignment[1].length(); i += 3 ) {
    
    aminoStr0 = alignment[0].substr( i, 3 );
    aminoStr1 = alignment[1].substr( i, 3 );
    aminoStr2 = alignment[2].substr( i, 3 );
    
    gapFlag = false;
    
    for ( int k = 0; k < 3; k++ ) {
      if ( aminoStr0[k] == '-' && aminoStr2[k] != '-') {
        data.insertion++;
	gapFlag = true;
      }
      
      if ( aminoStr0[k] != '-' && aminoStr2[k] == '-') {
        data.deletion++;
	gapFlag = true;
      }
    }
    
    if ( !gapFlag ) {
      bool synonymous = aminoCheck( aminoStr0, aminoStr2 ); 
    
      if ( synonymous ) data.synMutation++;
      else data.nonSynMutation++;
    }
    
    for ( int j = 0; j < 3; j++ ) {
      if ( aminoStr1[j] == '|' ) data.match++;
    }
    
    
  }
  
  data.mismatch = alignment[1].length() - ( data.match 
		  + data.insertion + data.deletion );
  
  return data;
}

