#ifndef SLENDERBODY_UTILS
#define SLENDERBODY_UTILS

#include <sctl.hpp>

namespace sctl {

template <class Real> class SlenderVolumeVis {
    static constexpr Integer COORD_DIM = 3;
    static constexpr Integer Nlevels = 15;
  public:

    SlenderVolumeVis() = default;

    SlenderVolumeVis(const Vector<Long>& cheb_order_, const Vector<Long>& fourier_order_, const Vector<Real>& coord_, const Vector<Real>& radius_, const Vector<Real>& orientation_ = Vector<Real>(), Integer cheb_upsample = 1, Integer fourier_upsample = 1) {
      const auto concat_vecs = [](Vector<Real>& v, const Vector<Vector<Real>>& vec_lst) {
        const Long N = vec_lst.Dim();
        Vector<Long> dsp(N+1); dsp[0] = 0;
        for (Long i = 0; i < N; i++) {
          dsp[i+1] = dsp[i] + vec_lst[i].Dim();
        }
        if (v.Dim() != dsp[N]) v.ReInit(dsp[N]);
        for (Long i = 0; i < N; i++) {
          Vector<Real> v_(vec_lst[i].Dim(), v.begin()+dsp[i], false);
          v_ = vec_lst[i];
        }
      };
      fourier_order = fourier_order_ * fourier_upsample;
      cheb_order = cheb_order_ * cheb_upsample;

      Vector<Vector<Real>> X;
      SlenderElemList<Real> elem_lst;
      Vector<Real> radius = radius_*0.99;
      for (Long i = 0; i < Nlevels; i++) {
        Vector<Real> X_, s_param, sin_theta, cos_theta;
        elem_lst.Init(cheb_order, fourier_order, coord_, radius, orientation_);
        for (Long j = 0; j < cheb_order.Dim(); j++) {
          s_param.ReInit(cheb_order[j]);
          sin_theta.ReInit(fourier_order[j]);
          cos_theta.ReInit(fourier_order[j]);
          for (Long k = 0; k < cheb_order[j]; k++) {
            s_param[k] = 0.5 - cos<Real>(const_pi<Real>()*k/(cheb_order[j]-1))*0.5;
          }
          for (Long k = 0; k < fourier_order[j]; k++) {
            sin_theta[k] = sin<Real>(2*const_pi<Real>()*k/fourier_order[j]);
            cos_theta[k] = cos<Real>(2*const_pi<Real>()*k/fourier_order[j]);
          }
          elem_lst.GetGeom(&X_, nullptr, nullptr, nullptr, nullptr, s_param, sin_theta, cos_theta, j);
          X.PushBack(X_);
        }
        radius *= 0.95;
      }
      concat_vecs(coord, X);
    }

    const Vector<Real>& GetCoord() const {
      return coord;
    }

    void GetVTUData(VTUData& vtu_data, const Vector<Real>& F) const {
      Long level_offset = 0;
      Vector<Long> elem_offset(cheb_order.Dim());
      for (Long i = 0; i < cheb_order.Dim(); i++) {
        elem_offset[i] = level_offset;
        level_offset += cheb_order[i] * fourier_order[i];
      }

      for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
      for (const auto& x :     F) vtu_data.value.PushBack((float)x);
      for (Long i = 0; i < cheb_order.Dim(); i++) {
        for (Long l = 0; l < Nlevels-1; l++) {
          for (Long j = 0; j < cheb_order[i]-1; j++) {
            for (Long k = 0; k < fourier_order[i]; k++) {
              auto idx = [&level_offset,&elem_offset,this](Integer elem_idx, Integer level, Integer cheb_idx, Integer fourier_idx) {
                return elem_offset[elem_idx] + level_offset*level + cheb_idx*fourier_order[elem_idx] + fourier_idx%fourier_order[elem_idx];
              };
              vtu_data.connect.PushBack(idx(i, l+0,j+0,k+0));
              vtu_data.connect.PushBack(idx(i, l+0,j+0,k+1));
              vtu_data.connect.PushBack(idx(i, l+0,j+1,k+1));
              vtu_data.connect.PushBack(idx(i, l+0,j+1,k+0));
              vtu_data.connect.PushBack(idx(i, l+1,j+0,k+0));
              vtu_data.connect.PushBack(idx(i, l+1,j+0,k+1));
              vtu_data.connect.PushBack(idx(i, l+1,j+1,k+1));
              vtu_data.connect.PushBack(idx(i, l+1,j+1,k+0));
              vtu_data.offset.PushBack(vtu_data.connect.Dim());;
              vtu_data.types.PushBack(12);
            }
          }
        }
      }
    }
    void WriteVTK(const std::string& fname, const Vector<Real>& F, const Comm& comm) const {
      VTUData vtu_data;
      GetVTUData(vtu_data, F);
      vtu_data.WriteVTK(fname, comm);
    }

  private:

    Vector<Long> cheb_order;
    Vector<Long> fourier_order;
    Vector<Real> coord;
};

}

#endif
