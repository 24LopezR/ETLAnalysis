#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVector3.h>
#include <iostream>
#include <string>
#include <TH1F.h>
#include <iomanip>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include <map>
#include <cmath>
#include <set>
#include <vector>
#include <regex>

#include <ROOT/RVec.hxx>
#include <set>

bool isGoodEvent(ROOT::VecOps::RVec<int> board) {
    if (board.size() != 4) return false;
    if (std::set(board.begin(),board.end()).size() != board.size()) return false;
    return true;
}

ROOT::VecOps::RVec<float> toaToPs(ROOT::VecOps::RVec<float> toa, ROOT::VecOps::RVec<int> acal) {
  ROOT::VecOps::RVec<float> toaps(toa.size());
  for (int i = 0; i < toa.size(); i++) 
    toaps[i] = 1000*(12.5-toa[i]*3.125/acal[i]);
  return toaps;
}

ROOT::VecOps::RVec<float> totToPs(ROOT::VecOps::RVec<float> tot, ROOT::VecOps::RVec<int> acal) {
  ROOT::VecOps::RVec<float> totps(tot.size());
  for (int i = 0; i < tot.size(); i++) 
    totps[i] = 1000.*((2.*tot[i]-floor(tot[i]/32.))*3.125/acal[i]);
  return totps;
}

bool passCalDiff(ROOT::VecOps::RVec<float> caldiff, float thresh=4.) {
  for (int i = 0; i < caldiff.size(); i++) {
    if (std::fabs(caldiff[i]) > thresh) return false;
  }
  return true;
}

bool isInPixel(ROOT::VecOps::RVec<int> board,
               ROOT::VecOps::RVec<int> col,
               ROOT::VecOps::RVec<int> row,
               int BOARD, int COL, int ROW) {
  for (int i = 0; i < board.size(); i++) {
    if (board[i] == BOARD) {
      if (col[i] == COL and row[i] == ROW) return true;
    }
  }
  return false;
}
