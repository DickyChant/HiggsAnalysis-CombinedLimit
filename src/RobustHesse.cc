#include "../interface/RobustHesse.h"

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>

#include "TH2F.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "RooWorkspace.h"
#include "TDecompBK.h"
#include "TMatrixDSymEigen.h"


RobustHesse::RobustHesse(RooAbsReal &nll, unsigned verbose) : nll_(&nll), verbosity_(verbose) {
  targetNllForStencils_ = 0.1;
  minNllForStencils_ = 0.095;
  maxNllForStencils_ = 0.105;
  doRescale_ = true;
  maxRemovalsFromHessian_ = 20;
  initialize();
}

void RobustHesse::initialize() {
  nll0_ = nll_->getVal();

  // Get a list of the floating RooRealVars
  std::unique_ptr<RooArgSet> allpars(nll_->getParameters(RooArgSet()));
  std::vector<Var> allVars;
  for (RooAbsArg *item : *allpars) {
    RooRealVar *rrv = dynamic_cast<RooRealVar*>(item);
    if (rrv && !rrv->isConstant()) {
      allVars.push_back(Var());
      allVars.back().v = rrv;
      allVars.back().nominal = rrv->getVal();
      // rrv->Print();
    }
  }
  std::cout << ">> Found " << allVars.size() << " floating parameters\n";

  ReplaceVars(allVars);
}

double RobustHesse::deltaNLL() {
  return nll_->getVal() - nll0_;
}

double RobustHesse::deltaNLL(std::vector<unsigned> const& indices, std::vector<double> const& vals) {
  if (indices.size() != vals.size()) {
    std::cout << ">> Error in deltaNLL: number of indices != number of vals\n";
    return 0.;
  }
  nllEvals_++;
  if (indices.size() == 1) {
    auto it = nllcache_.find(std::make_pair(indices[0], vals[0]));
    if (it == nllcache_.end()) {
      for (unsigned i = 0; i < indices.size(); ++i) {
        cVars_[indices[i]].v->setVal(vals[i]);
      }
      it = nllcache_.emplace(std::make_pair(indices[0], vals[0]), deltaNLL()).first;
      for (unsigned i = 0; i < indices.size(); ++i) {
        cVars_[indices[i]].v->setVal(cVars_[indices[i]].nominal);
      }
    } else {
      nllEvalsCached_++;
    }
    return it->second;
  }
  for (unsigned i = 0; i < indices.size(); ++i) {
    cVars_[indices[i]].v->setVal(vals[i]);
  }
  double result = deltaNLL();

  for (unsigned i = 0; i < indices.size(); ++i) {
    cVars_[indices[i]].v->setVal(cVars_[indices[i]].nominal);
  }
  return result;
}


int RobustHesse::setParameterStencil(unsigned i) {
  double x = cVars_[i].nominal;
  RooRealVar *rrv = cVars_[i].v;

  double valLo = x - rrv->getError();
  double valHi = x + rrv->getError();

  // Am I near a boundary?
  double boundLo = rrv->getMin();
  double boundHi = rrv->getMax();

  // If we skip initial fit (e.g. to compute from asimov point)
  // need a guess for the initial step size
  if ( rrv->getError() == 0.0 ) {
    valLo = x + 1e-3*(boundLo - x);
    valHi = x + 1e-3*(boundHi - x);
  }

  bool closeToLo = valLo < boundLo;
  bool closeToHi = valHi > boundHi;

  if (closeToLo) {
    if (verbosity_ > 0) {
      std::cout << ">> Parameter " << rrv->GetName() << " is too close to the lower bound: \n";
      rrv->Print();
    }
    valLo = boundLo + 1E-7;
  }
  if (closeToHi) {
    if (verbosity_ > 0) {
      std::cout << ">> Parameter " << rrv->GetName() << " is too close to the upper bound: \n";
      rrv->Print();
    }
    valHi = boundHi + 1E-7;
  }

  double dNllLo = deltaNLL({i}, {valLo});
  double dNllHi = deltaNLL({i}, {valHi});

  // What sort of problems can we have here?
  //  a) dNll can be negative => we weren't quite at the minimum,
  //     but maybe we did just take a very small step
  //  b) dNLL is too large or too small. Should perform a search to get in the target range
  //  c) dNLL could be zero, or zero-ish. Might be a non-active parameter in the multipdf
  // Solution to cases (a) and (b) is the same: perform a search in the allowed range to see if
  // we can find a better val. If we don't on one side only, then we can use an asymmetric stencil
  // if we fail on both sides we have a problem.
  // Note that a very large dNll could just mean the parameter has a very strong anti-correlation
  // with another.

  // std::cout << rrv->GetName() << "\t" << dNllLo << "\t" << dNllHi << "\n";

  bool notViableLo = false;
  bool notViableHi = false;

  if (dNllLo < minNllForStencils_) {
    if (closeToLo) {
      notViableLo = true;
    } else {
      auto bound_result = findBound(i, x, valLo - x, 2., 2., targetNllForStencils_, boundLo, 100);
      valLo = bound_result.second;
      dNllLo = deltaNLL({i}, {valLo});
      if (dNllLo < minNllForStencils_) {
        if (verbosity_ > 0) std::cout << ">> dNllLo is still too small after [findBound]\n";
        notViableLo = true;
      } else {
        if (verbosity_ > 0) std::cout << ">> dNllLo is now " << dNllLo << "\n";
      }
    }
    if (notViableLo && verbosity_ > 0) std::cout << ">> dNllLo is not viable\n";
  }

  if (dNllHi < minNllForStencils_) {
    if (closeToHi) {
      notViableHi = true;
    } else {
      auto bound_result = findBound(i, x, valHi - x, 2., 2., targetNllForStencils_, boundHi, 100);
      valHi = bound_result.second;
      dNllHi = deltaNLL({i}, {valHi});
      if (dNllHi < minNllForStencils_) {
        if (verbosity_ > 0) std::cout << ">> dNllHi is still too small after [findBound]\n";
        notViableHi = true;
      } else {
        if (verbosity_ > 0) std::cout << ">> dNllHi is now " << dNllHi << "\n";
      }
    }
    if (notViableHi && verbosity_ > 0) std::cout << ">> dNllHi is not viable\n";
  }

  if (dNllLo > maxNllForStencils_) {
    valLo = improveWithBisect(i, x, valLo, targetNllForStencils_, minNllForStencils_,  maxNllForStencils_, 20);
    dNllLo = deltaNLL({i}, {valLo});
  }
  if (dNllHi > maxNllForStencils_) {
    valHi = improveWithBisect(i, x, valHi, targetNllForStencils_, minNllForStencils_,  maxNllForStencils_, 20);
    dNllHi = deltaNLL({i}, {valHi});
  }



  // std::cout << "[setParameterStencil] " << rrv->GetName() << "\t" << x << "\t" << valLo << " (" << dNllLo << ")\t" << valHi << " (" << dNllHi << ")\n";
  if (verbosity_ > 0) printf("%-80s %-10.3f | %-10.3f %-10.3f %-4.i | %-10.3f (%-10.3f) %-4.i \n",
    rrv->GetName(), x, valLo, dNllLo, !notViableLo, valHi, dNllHi, !notViableHi);

  cVars_[i].rescale = 1.0;
  if (!notViableLo && !notViableHi) {
    if (doRescale_) {
      cVars_[i].rescale = (valHi - valLo) / 2.;
    }
    cVars_[i].stencil = {(valLo - x) / cVars_[i].rescale, 0., (valHi - x) / cVars_[i].rescale};
  } else if (notViableLo && !notViableHi) {
    if (doRescale_) {
      cVars_[i ].rescale = (valHi - x);
    }
    cVars_[i].stencil = {0.,  (valHi - x) / (2. * cVars_[i].rescale) ,(valHi - x) / cVars_[i].rescale};
  } else if (!notViableLo && notViableHi) {
    if (doRescale_) {
      cVars_[i].rescale = (x - valLo);
    }
    cVars_[i].stencil = {(valLo - x) / cVars_[i].rescale, (valLo - x) / (2. * cVars_[i].rescale), 0.};
  } else {
    std::cout << ">> Parameter " << rrv->GetName() << " does not have a valid stencil and will be dropped\n";
  }

  return 0;
}

std::pair<int, double> RobustHesse::findBound(unsigned i, double x, double initialDelta, double initialMult, double scaleMult, double threshold, double hardBound, unsigned maxIters) {
  if (verbosity_ > 0) std::cout << ">> [findBound] for parameter " << i << ", best-fit = " << x << ", initialDelta = " << initialDelta << "\n";
  double mult = initialMult;
  double trial = x + initialDelta;
  for (unsigned j = 0; j < maxIters; ++j) {
    bool trial_hit_bound = false;
    trial = x + initialDelta * mult;
    if (hardBound > x && trial > hardBound) {
      if (verbosity_ > 0) std::cout << "Step would go past bound, reducing\n";
      trial = hardBound - 1E-7;
      trial_hit_bound = true;
    }
    if (hardBound < x && trial < hardBound) {
      if (verbosity_ > 0) std::cout << "Step would go past bound, reducing\n";
      trial = hardBound + 1E-7;
      trial_hit_bound = true;
    }

    double dNll = deltaNLL({i}, {trial});
    if (verbosity_ > 0) std::cout << ">> dNLL at " << trial << " (x" << mult << " of original delta) = " << dNll << "\n";
    if ((mult > 1. && dNll > threshold) || (mult <= 1. && dNll < threshold)) {
      if (verbosity_ > 0) std::cout << ">> This is past the threshold, returning\n";
      return std::make_pair(0, trial);
    } else if (trial_hit_bound) {
      return std::make_pair(1, trial);
    }
    mult *= scaleMult;
  }
  return std::make_pair(2, trial);
}


double RobustHesse::improveWithBisect(unsigned i, double x, double max, double target, double targetLo, double targetHi, unsigned maxIters) {
  double cmin = x;
  double cmax = max;
  double trial = max;
  for (unsigned j = 0; j < maxIters; ++j) {
    trial = (cmin + cmax) / 2.;
    double dNll = deltaNLL({i}, {trial});
    if (verbosity_ > 0) std::cout << "[Bisect] " << j << " cmin = " << cmin << ", cmax = " << cmax << ", trial = " << trial << " dNll = " << dNll << "\n";
    if (dNll > targetLo && dNll < targetHi) {
      break;
    }
    if (dNll < target) {
      cmin = trial;
    } else {
      cmax = trial;
    }
  }
  return trial;
}

std::vector<double> RobustHesse::getFdCoeffs(unsigned n, std::vector<double> const& stencil) {
  std::vector<double> result;
  unsigned N = stencil.size();
  TMatrixD smatrix(N, N);
  TVectorD dvec(N);
  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < N; ++j) {
      smatrix(i, j) = std::pow(stencil[j], double(i));
    }
    dvec[i] = double((i == n) ? factorial(n) : 0);
  }
  TMatrixD smatrix2 = smatrix;
  smatrix2.Invert();
  TVectorD res = (smatrix2 * dvec);
  for (unsigned i = 0; i < N; ++i) {
    result.push_back(res[i]);
  }
  return result;
}


void RobustHesse::RemoveFromHessian(std::vector<unsigned> const& ids) {

  std::unique_ptr<TMatrixDSym> hessian2ptr(new TMatrixDSym(cVars_.size() - ids.size()));
  TMatrixDSym & hessian2 = *(hessian2ptr.get());

  unsigned inew = 0;
  unsigned jnew = 0;
  std::vector<bool> ikeep(cVars_.size());

  for (unsigned i = 0; i < cVars_.size(); ++i) {
    ikeep[i] = (std::find(ids.begin(), ids.end(), i) == ids.end());
  }

  std::vector<Var> keptVars;
  for (unsigned i = 0; i < cVars_.size(); ++i) {
    if (ikeep[i]) {
      keptVars.push_back(cVars_[i]);
    } else {
      removedFromHessianVars_.push_back(cVars_[i]);
      std::cout << ">> Dropping " << cVars_[i].v->GetName() << " from the hessian\n";
    }
    jnew = inew;
    if (ikeep[i]) {
      for (unsigned j = i; j < cVars_.size(); ++j) {
        if (ikeep[j]) {
          hessian2[inew][jnew] = (*hessian_)[i][j];
          hessian2[jnew][inew] = (*hessian_)[j][i];
          ++jnew;
        }
      }
      ++inew;
    }
  }

  ReplaceVars(keptVars);
  hessian_ = std::move(hessian2ptr);
}



int RobustHesse::hesse(unsigned NsplitHesse, unsigned iHesse) {

  if (NsplitHesse < 1 || iHesse < 1 || iHesse > NsplitHesse) {
    std::cerr << ">> Error in RobustHesse::hesse: invalid split configuration iHesse = " << iHesse << ", NsplitHesse = " << NsplitHesse << " (need 1 <= iHesse <= NsplitHesse)\n";
    return 1;
  }
  const bool partial = (NsplitHesse > 1) && mergeFile_.empty() && loadFile_.empty();
  if (partial && saveFile_.empty()) {
    std::cerr << ">> Error in RobustHesse::hesse: a split Hessian job (NsplitHesse > 1) needs a save file (--robustHesseSave) for its part\n";
    return 1;
  }

  // Step 1: try and set parameter stencils at the target NLL values
  for (unsigned i = 0; i < cVars_.size(); ++i) {
    setParameterStencil(i);
  }

  std::vector<Var> validStencilVars;
  for (unsigned i = 0; i < cVars_.size(); ++i) {
    if (cVars_[i].stencil.size() == 3) {
      validStencilVars.push_back(cVars_[i]);
    } else {
      invalidStencilVars_.push_back(cVars_[i]);
    }
  }
  ReplaceVars(validStencilVars);

  for (unsigned i = 0; i < cVars_.size(); ++i) {
    if (cVars_[i].stencil.size() == 3) {
      cVars_[i].d1coeffs = getFdCoeffs(1, cVars_[i].stencil);
      cVars_[i].d2coeffs = getFdCoeffs(2, cVars_[i].stencil);
      if (verbosity_ > 0) printf("%-80s %-10.3f %-10.3f %-10.3f | %-10.3f %-10.3f %-10.3f | %-10.3f %-10.3f %-10.3f\n",
        cVars_[i].v->GetName(), cVars_[i].stencil[0], cVars_[i].stencil[1], cVars_[i].stencil[2],
        cVars_[i].d1coeffs[0], cVars_[i].d1coeffs[1], cVars_[i].d1coeffs[2],
        cVars_[i].d2coeffs[0], cVars_[i].d2coeffs[1], cVars_[i].d2coeffs[2]
        );
    }
  }

  // Step 2: Calculate and populate hessian
  const unsigned nvars = cVars_.size();
  hessian_ = std::unique_ptr<TMatrixDSym>(new TMatrixDSym(nvars));
  unsigned ntotal = (((nvars * nvars) - nvars) / 2) + nvars;
  unsigned idx = 0;

  if (!mergeFile_.empty()) {
    if (!mergePartialHessian()) return 1;
  } else if (loadFile_ != "") {
    TFile fin(loadFile_.c_str());
    if (fin.IsZombie() || !fin.IsOpen()) {
      std::cerr << ">> Error in RobustHesse::hesse: cannot open Hessian file " << loadFile_ << "\n";
      return 1;
    }
    TMatrixDSym *loaded = dynamic_cast<TMatrixDSym*>(fin.Get("hessian"));
    if (!loaded) {
      std::cerr << ">> Error in RobustHesse::hesse: no TMatrixDSym \"hessian\" in " << loadFile_ << "\n";
      return 1;
    }
    if (loaded->GetNrows() != int(nvars)) {
      std::cerr << ">> Error in RobustHesse::hesse: Hessian in " << loadFile_ << " has dimension " << loaded->GetNrows()
                << " but " << nvars << " parameters have valid stencils in this job\n";
      return 1;
    }
    *hessian_ = *loaded;
  } else {
    // Range [idxMin, idxMax) of the upper-triangle cells (row-major, j >= i) computed by this job
    unsigned idxMin = 0;
    unsigned idxMax = ntotal;
    if (NsplitHesse > 1) {
      if (NsplitHesse > ntotal) {
        std::cerr << ">> Error in RobustHesse::hesse: NsplitHesse (" << NsplitHesse << ") must not exceed the number of independent Hessian cells ("
                  << ntotal << " for " << nvars << " parameters)\n";
        return 1;
      }
      unsigned cellsPerJob = ntotal / NsplitHesse;
      idxMin = cellsPerJob * (iHesse - 1);
      idxMax = cellsPerJob * iHesse;
      if (iHesse == NsplitHesse) idxMax = ntotal;  // last job takes the remainder
      std::cout << ">> Split Hessian: job " << iHesse << "/" << NsplitHesse << " computes cells [" << idxMin << ", " << idxMax << ") of " << ntotal << "\n";
    }
    for (unsigned i = 0; i < nvars; ++i) {
      for (unsigned j = i; j < nvars; ++j) {
        if (idx < idxMin || idx >= idxMax) {
          ++idx;
          continue;
        }
        if (idx % 100 == 0) {
          if (verbosity_ > 0) std::cout << " - Done " << idx << "/" << ntotal << " terms (" << nllEvals_ << " evals, of which " << nllEvalsCached_ << " cached)\n";
        }
        double term = 0.;
        if (i == j) {
          for (unsigned k = 0; k < cVars_[i].stencil.size(); ++k) {
            if (cVars_[i].stencil[k] != 0.) {
              term += deltaNLL({i}, {cVars_[i].nominal + cVars_[i].rescale * cVars_[i].stencil[k]}) * cVars_[i].d2coeffs[k];
            }
          }
        } else {
          for (unsigned k = 0; k < cVars_[i].stencil.size(); ++k) {
            double c1 = cVars_[i].d1coeffs[k];
            double c2 = 0.;
            for (unsigned l = 0; l < cVars_[j].stencil.size();++l) {
              if (cVars_[i].stencil[k] == 0. && cVars_[j].stencil[l] == 0.) {
                continue;
              } else if (cVars_[i].stencil[k] == 0.) {
                c2 += deltaNLL({j}, {cVars_[j].nominal + cVars_[j].stencil[l] * cVars_[j].rescale}) * cVars_[j].d1coeffs[l];
              } else if (cVars_[j].stencil[l] == 0.) {
                c2 += deltaNLL({i}, {cVars_[i].nominal + cVars_[i].stencil[k] * cVars_[i].rescale}) * cVars_[j].d1coeffs[l];
              } else {
                c2 += deltaNLL({i, j}, {cVars_[i].nominal + cVars_[i].stencil[k] * cVars_[i].rescale, cVars_[j].nominal + cVars_[j].stencil[l] * cVars_[j].rescale}) * cVars_[j].d1coeffs[l];
              }
            }
            term += (c2 * c1);
          }
        }
        (*hessian_)[i][j] = term;
        (*hessian_)[j][i] = term;
        ++idx;
      }
    }
    if (NsplitHesse > 1) {
      std::string partFile = PartialHessianFileName(saveFile_, iHesse, NsplitHesse);
      if (!savePartialHessian(partFile, NsplitHesse, iHesse, idxMin, idxMax, ntotal)) return 1;
      std::cout << ">> Partial Hessian (job " << iHesse << "/" << NsplitHesse << ") saved to " << partFile << "\n";
      return 0;
    }
    if (saveFile_ != "") {
      TFile fout(saveFile_.c_str(), "RECREATE");
      gDirectory->WriteObject(hessian_.get(), "hessian");
    }
  }

  bool print_only_negative = true;
  std::vector<std::string> removed_pars;
  for (unsigned ai = 0; ai < maxRemovalsFromHessian_; ++ai) {
    std::cout << ">> Verifying eigenvalues are all positive (iteration " << ai << ")\n";
    bool have_negative_eigenvalues = false;

    TMatrixDSymEigen eigen(*hessian_);
    auto eigenvals = eigen.GetEigenValues();
    auto eigenvecs = eigen.GetEigenVectors();

    std::vector<unsigned> to_remove;

    for (int ei = 0; ei < eigenvals.GetNrows(); ++ei) {
      if (eigenvals[ei] < 0.) have_negative_eigenvalues = true;

      std::vector<std::pair<unsigned, double>> eigenvec(eigenvecs.GetNrows());
      for (int ej = 0; ej < eigenvecs.GetNrows(); ++ej) {
        eigenvec[ej] = std::make_pair(unsigned(ej), eigenvecs[ej][ei]);
      }
      if ((!print_only_negative) || (print_only_negative && eigenvals[ei] < 0.)) {
        std::cout << "Eigenvalue " << ei << " : " << eigenvals[ei] << "\n";
        std::sort(eigenvec.begin(), eigenvec.end(), [](std::pair<unsigned, double> const& p1, std::pair<unsigned, double> const& p2) {
          return std::fabs(p1.second) > std::fabs(p2.second);
        });
        unsigned n_show = 5;
        for (unsigned es = 0; es < TMath::Min(n_show, unsigned(eigenvec.size())); ++es) {
          std::cout << " - " << eigenvec[es].second << "\t" << cVars_[eigenvec[es].first].v->GetName() << "\n";
        }

        
        // Need to ensure this only gets triggered for negative eigenvalues in case print_only_negative is set to false
        if (eigenvals[ei] < 0.){
          // Assume sorted by largest to smallest eigenvalues,
          // meaning most negative is this one
          if (ei == (eigenvals.GetNrows() - 1)) {
            for (unsigned ej = 0; ej < eigenvec.size(); ++ ej) {
              // Check if we are allowed to remove this parameter
              std::string parname = cVars_[eigenvec[ej].first].v->GetName();
              // If it is not protected OR if it is the last parameter then we can remove it
              if (!proctected_.count(parname) || ej == (eigenvec.size() - 1)) {
                to_remove.push_back(eigenvec[ej].first);
                removed_pars.push_back(parname);
                break;
              }
            }
          }
        }
      }
    }
    if (!have_negative_eigenvalues) {
      std::cout << ">> All eigenvalues are positive\n";
      break;
    } else {
      std::cout << ">> Attempt " << ai << " failed, remove the variable with largest eigenvector component...\n";
      RemoveFromHessian(to_remove);
    }
  }

  covariance_ = std::unique_ptr<TMatrixDSym>(new TMatrixDSym(*hessian_));
  std::cout << ">> Inverting Hessian...\n";
  covariance_->Invert();
  std::cout << ">> Matrix inverted\n";

  if (invalidStencilVars_.size() > 0) {
    std::cout << ">> Note: the following parameters were dropped as their ranges are either too narrow or the parameter has no effect on the NLL:\n";
    for (auto const& par : invalidStencilVars_) {
      std::cout << "     - " << par.v->GetName() << "\n";
    }
  }
  if (removed_pars.size() > 0) {
    std::cout << ">> Note: the following parameters were dropped to make the covariance matrix pos-def:\n";
    for (auto const& par : removed_pars) {
      std::cout << "     - " << par << "\n";
    }
  }

  // TDecompBK decomp(hessian);
  // std::cout << " - Inverting Hessian...\n";
  // TMatrixDSym covariance = decomp.Invert();
  // covariance.Invert();


  if (verbosity_ > 0) {
    for (unsigned i = 0; i < cVars_.size(); ++i) {
      std::cout << cVars_[i].v->GetName() << "\t" << (std::sqrt((*covariance_)[i][i]) * cVars_[i].rescale) << "\n";
      for (unsigned j = i+1; j < cVars_.size(); ++j) {
        double corr = (*covariance_)[i][j] / (std::sqrt((*covariance_)[i][i]) * std::sqrt((*covariance_)[j][j]));
        if (std::abs(corr) > 0.8) {
          std::cout << " - correlation with " << cVars_[j].v->GetName() << " = " << corr << "\n";
        }
      }
    }
  }
  return 0;
}

void RobustHesse::WriteOutputFile(std::string const& outputFileName) const {
  TFile fout(outputFileName.c_str(), "RECREATE");
  TH2F h_cov("covariance", "covariance", cVars_.size(), 0, cVars_.size(), cVars_.size(), 0,
             cVars_.size());
  TH2F h_cor("correlation", "correlation", cVars_.size(), 0, cVars_.size(), cVars_.size(), 0,
             cVars_.size());
  RooArgList arglist("floatParsFinal");
  for (unsigned i = 0; i < cVars_.size(); ++i) {
    arglist.addClone(*cVars_[i].v);
    RooRealVar *rrv = dynamic_cast<RooRealVar*>(arglist.at(i));
    rrv->setError(std::sqrt((*covariance_)[i][i]) * cVars_[i].rescale);
    h_cov.GetXaxis()->SetBinLabel(i + 1, cVars_[i].v->GetName());
    h_cov.GetYaxis()->SetBinLabel(i + 1, cVars_[i].v->GetName());
    h_cor.GetXaxis()->SetBinLabel(i + 1, cVars_[i].v->GetName());
    h_cor.GetYaxis()->SetBinLabel(i + 1, cVars_[i].v->GetName());
    for (unsigned j = i; j < cVars_.size(); ++j) {
      h_cov.SetBinContent(i + 1, j + 1,
                          (*covariance_)[i][j] * cVars_[i].rescale * cVars_[j].rescale);
      h_cov.SetBinContent(j + 1, i + 1,
                          (*covariance_)[i][j] * cVars_[i].rescale * cVars_[j].rescale);
      h_cor.SetBinContent(i + 1, j + 1,
                          (*covariance_)[i][j] / (std::sqrt((*covariance_)[i][i]) * std::sqrt((*covariance_)[j][j])));
      h_cor.SetBinContent(j + 1, i + 1,
                          (*covariance_)[i][j] / (std::sqrt((*covariance_)[i][i]) * std::sqrt((*covariance_)[j][j])));
    }
  }
  gDirectory->WriteObject(hessian_.get(), "hessian");
  gDirectory->WriteObject(covariance_.get(), "covariance");
  gDirectory->WriteObject(&h_cov, "h_covariance");
  gDirectory->WriteObject(&h_cor, "h_correlation");
  gDirectory->WriteObject(arglist.snapshot(), "floatParsFinal");
  fout.Close();
}

RooFitResult* RobustHesse::GetRooFitResult(RooFitResult const* current) const {
  RooFitResultBuilder* rfr = nullptr;
  if (current) {
    rfr = new RooFitResultBuilder(*current);
  } else {
    rfr = new RooFitResultBuilder();
  }

  auto scaled_cov = std::unique_ptr<TMatrixDSym>(new TMatrixDSym(cVars_.size() + removedFromHessianVars_.size()));
  RooArgList arglist("floatParsFinal");

  for (unsigned i = 0; i < cVars_.size(); ++i) {
    RooRealVar newVar(cVars_[i].v->GetName(), "", cVars_[i].v->getVal(), cVars_[i].v->getMin(), cVars_[i].v->getMax());
    arglist.addClone(newVar);
    RooRealVar* rrv = dynamic_cast<RooRealVar*>(arglist.at(i));
    rrv->setError(std::sqrt((*covariance_)[i][i]) * cVars_[i].rescale);
    for (unsigned j = i; j < cVars_.size(); ++j) {
      (*scaled_cov)[i][j] = (*covariance_)[i][j] * cVars_[i].rescale * cVars_[j].rescale;
      (*scaled_cov)[j][i] = (*covariance_)[i][j] * cVars_[i].rescale * cVars_[j].rescale;
    }
  }
  for (unsigned i = 0; i < removedFromHessianVars_.size(); ++i) {
    RooRealVar newVar(
        removedFromHessianVars_[i].v->GetName(), "", removedFromHessianVars_[i].v->getVal(),
        removedFromHessianVars_[i].v->getMin(), removedFromHessianVars_[i].v->getMax());
    arglist.addClone(newVar);
    RooRealVar* rrv = dynamic_cast<RooRealVar*>(arglist.at(cVars_.size() + i));
    rrv->setError(0.);
    // Set some nominally small value
    (*scaled_cov)[cVars_.size() + i][cVars_.size() + i] = (1E-100);

  }
  rfr->setFinalParList(arglist);
  rfr->setCovarianceMatrix(*scaled_cov);
  return rfr;
}

void RobustHesse::SaveHessianToFile(std::string const& filename) {
  saveFile_ = filename;
}
void RobustHesse::LoadHessianFromFile(std::string const& filename) {
  loadFile_ = filename;
}

std::string RobustHesse::PartialHessianFileName(std::string const& filename, unsigned iHesse, unsigned NsplitHesse) {
  std::string label = ".part" + std::to_string(iHesse) + "of" + std::to_string(NsplitHesse);
  const std::string ext = ".root";
  if (filename.size() >= ext.size() && filename.compare(filename.size() - ext.size(), ext.size(), ext) == 0) {
    return filename.substr(0, filename.size() - ext.size()) + label + ext;
  }
  return filename + label + ext;
}

bool RobustHesse::savePartialHessian(std::string const& filename, unsigned NsplitHesse, unsigned iHesse, unsigned idxMin, unsigned idxMax, unsigned ntotal) const {
  TFile fout(filename.c_str(), "RECREATE");
  if (fout.IsZombie() || !fout.IsOpen()) {
    std::cerr << ">> Error in RobustHesse: cannot create partial Hessian file " << filename << "\n";
    return false;
  }
  // Cells computed by this job, both (i,j) and (j,i) (diagonal once)
  unsigned row = 0, col = 0;
  double hessianval = 0.;
  fout.cd();
  // Trees are owned by (and deleted with) the output file
  TTree *tcells = new TTree("partialhessian", "partialhessian");
  tcells->Branch("row", &row, "row/i");
  tcells->Branch("col", &col, "col/i");
  tcells->Branch("hessianval", &hessianval, "hessianval/D");
  unsigned idx = 0;
  for (unsigned i = 0; i < cVars_.size(); ++i) {
    for (unsigned j = i; j < cVars_.size(); ++j, ++idx) {
      if (idx < idxMin || idx >= idxMax) continue;
      row = i; col = j; hessianval = (*hessian_)[i][j];
      tcells->Fill();
      if (i != j) {
        row = j; col = i;
        tcells->Fill();
      }
    }
  }
  // Parameters (in Hessian order) with their nominal values and stencil rescale factors,
  // used at merge time to check that all parts are compatible with each other
  std::string name;
  double nominal = 0., rescale = 0.;
  TTree *tpars = new TTree("hessianpars", "hessianpars");
  tpars->Branch("name", &name);
  tpars->Branch("nominal", &nominal, "nominal/D");
  tpars->Branch("rescale", &rescale, "rescale/D");
  for (auto const& var : cVars_) {
    name = var.v->GetName();
    nominal = var.nominal;
    rescale = var.rescale;
    tpars->Fill();
  }
  TVectorD info(6);
  info[0] = NsplitHesse;
  info[1] = iHesse;
  info[2] = idxMin;
  info[3] = idxMax;
  info[4] = ntotal;
  info[5] = cVars_.size();
  fout.cd();
  tcells->Write();
  tpars->Write();
  info.Write("splitinfo");
  fout.Close();
  return true;
}

bool RobustHesse::MergePartialHessianFromFile(std::string const& filename, unsigned NsplitHesse) {
  if (NsplitHesse < 2) {
    std::cerr << ">> Error in RobustHesse: merging partial Hessians requires NsplitHesse > 1 (got " << NsplitHesse << ")\n";
    return false;
  }
  for (unsigned ipart = 1; ipart <= NsplitHesse; ++ipart) {
    std::string partFile = PartialHessianFileName(filename, ipart, NsplitHesse);
    TFile fin(partFile.c_str(), "READ");
    if (fin.IsZombie() || !fin.IsOpen() || !fin.Get("partialhessian")) {
      std::cerr << ">> Error in RobustHesse: partial Hessian file " << partFile << " is missing or invalid\n";
      return false;
    }
  }
  mergeFile_ = filename;
  mergeNsplit_ = NsplitHesse;
  return true;
}

bool RobustHesse::mergePartialHessian() {
  const unsigned n = cVars_.size();
  const unsigned N = mergeNsplit_;
  std::vector<char> seen(size_t(n) * n, 0);
  TMatrixDSym merged(n);
  double maxNominalDiff = 0., maxRescaleDiff = 0.;
  std::string maxNominalPar, maxRescalePar;
  size_t ncells = 0;

  for (unsigned ipart = 1; ipart <= N; ++ipart) {
    std::string partFile = PartialHessianFileName(mergeFile_, ipart, N);
    TFile fin(partFile.c_str(), "READ");
    if (fin.IsZombie() || !fin.IsOpen()) {
      std::cerr << ">> Error in RobustHesse: cannot open partial Hessian file " << partFile << "\n";
      return false;
    }
    TTree *tcells = dynamic_cast<TTree*>(fin.Get("partialhessian"));
    TTree *tpars = dynamic_cast<TTree*>(fin.Get("hessianpars"));
    TVectorD *info = dynamic_cast<TVectorD*>(fin.Get("splitinfo"));
    if (!tcells || !tpars || !info || info->GetNrows() < 6) {
      std::cerr << ">> Error in RobustHesse: " << partFile << " does not contain a valid partial Hessian (partialhessian, hessianpars, splitinfo)\n";
      return false;
    }
    if (unsigned((*info)[0]) != N || unsigned((*info)[1]) != ipart) {
      std::cerr << ">> Error in RobustHesse: " << partFile << " was produced as part " << unsigned((*info)[1]) << " of " << unsigned((*info)[0])
                << ", expected part " << ipart << " of " << N << "\n";
      return false;
    }
    // Parameter list must match the current job exactly (same names, same order)
    if (unsigned((*info)[5]) != n || unsigned(tpars->GetEntries()) != n) {
      std::cerr << ">> Error in RobustHesse: " << partFile << " has a Hessian of dimension " << unsigned((*info)[5])
                << " but " << n << " parameters have valid stencils in this job\n";
      return false;
    }
    std::string *pname = nullptr;
    double nominal = 0., rescale = 0.;
    tpars->SetBranchAddress("name", &pname);
    tpars->SetBranchAddress("nominal", &nominal);
    tpars->SetBranchAddress("rescale", &rescale);
    for (unsigned k = 0; k < n; ++k) {
      tpars->GetEntry(k);
      if (!pname || *pname != cVars_[k].v->GetName()) {
        std::cerr << ">> Error in RobustHesse: parameter " << k << " in " << partFile << " is " << (pname ? *pname : std::string("?"))
                  << " but " << cVars_[k].v->GetName() << " in this job (different workspace or parameter ordering)\n";
        return false;
      }
      double dnom = std::fabs(nominal - cVars_[k].nominal) / std::max(1., std::fabs(cVars_[k].nominal));
      double dres = std::fabs(rescale - cVars_[k].rescale) / std::max(1e-300, std::fabs(cVars_[k].rescale));
      if (dnom > maxNominalDiff) { maxNominalDiff = dnom; maxNominalPar = cVars_[k].v->GetName(); }
      if (dres > maxRescaleDiff) { maxRescaleDiff = dres; maxRescalePar = cVars_[k].v->GetName(); }
    }
    tpars->ResetBranchAddresses();

    unsigned row = 0, col = 0;
    double hessianval = 0.;
    tcells->SetBranchAddress("row", &row);
    tcells->SetBranchAddress("col", &col);
    tcells->SetBranchAddress("hessianval", &hessianval);
    for (Long64_t ientry = 0; ientry < tcells->GetEntries(); ++ientry) {
      tcells->GetEntry(ientry);
      if (row >= n || col >= n) {
        std::cerr << ">> Error in RobustHesse: cell (" << row << ", " << col << ") in " << partFile << " is outside the " << n << "x" << n << " Hessian\n";
        return false;
      }
      char &s = seen[size_t(row) * n + col];
      if (s) {
        std::cerr << ">> Error in RobustHesse: cell (" << row << ", " << col << ") appears more than once (duplicate in " << partFile << ")\n";
        return false;
      }
      s = 1;
      merged[row][col] = hessianval;
      ++ncells;
    }
    tcells->ResetBranchAddresses();
  }

  if (ncells != size_t(n) * n) {
    size_t nmissing = 0;
    for (unsigned i = 0; i < n; ++i) {
      for (unsigned j = 0; j < n; ++j) {
        if (!seen[size_t(i) * n + j]) {
          if (nmissing < 10) std::cerr << ">> Error in RobustHesse: Hessian cell (" << i << ", " << j << ") missing from the partial files\n";
          ++nmissing;
        }
      }
    }
    std::cerr << ">> Error in RobustHesse: " << nmissing << " of " << size_t(n) * n << " Hessian cells missing after merging " << N << " parts\n";
    return false;
  }
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = i + 1; j < n; ++j) {
      if (merged[i][j] != merged[j][i]) {
        std::cerr << ">> Error in RobustHesse: merged Hessian is not symmetric at (" << i << ", " << j << ")\n";
        return false;
      }
    }
  }
  const double tol = 1e-6;
  if (maxNominalDiff > tol || maxRescaleDiff > tol) {
    std::cout << ">> Warning in RobustHesse: the partial Hessians were not computed at the same point as this job"
              << " (max rel. diff of nominal value " << maxNominalDiff << " for " << maxNominalPar
              << ", of stencil rescale " << maxRescaleDiff << " for " << maxRescalePar << ")."
              << " All split jobs should start from the same minimum, e.g. --snapshotName MultiDimFit --skipInitialFit\n";
  }
  std::cout << ">> Merged " << N << " partial Hessians (" << ncells << " cells, " << n << " parameters) into " << mergeFile_ << "\n";

  *hessian_ = merged;
  TFile fout(mergeFile_.c_str(), "RECREATE");
  if (fout.IsZombie() || !fout.IsOpen()) {
    std::cerr << ">> Error in RobustHesse: cannot write merged Hessian to " << mergeFile_ << "\n";
    return false;
  }
  gDirectory->WriteObject(hessian_.get(), "hessian");
  fout.Close();
  return true;
}

void RobustHesse::ProtectArgSet(RooArgSet const& set) {
  std::vector<std::string> names;
  for (RooAbsArg *item : set) {
    RooRealVar *rrv = dynamic_cast<RooRealVar*>(item);
    if (rrv && !rrv->isConstant()) {
      names.push_back(rrv->GetName());
    }
  }
  ProtectVars(names);
}

void RobustHesse::ProtectVars(std::vector<std::string> const& names) {
  std::cout << ">> The following parameters are protected in the event of negative eigenvalues:";
  for (auto const& name : names) {
    std::cout << " " << name;
    proctected_.insert(name);
  }
  std::cout << "\n";
}
