// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Mi Jiang
//
// Three-band model added into s orbital
// 260106Peng: 3D d-p-s model

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::dps_3D_model<PointGroup>>> {
public:
  ModelParameters() {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t_pd() const {
    return t_pd_;
  }
  void set_t_pd(const double t_pd) {
    t_pd_ = t_pd;
  }

  double get_t_ps() const {                                                        
    return t_ps_; 
  }
  void set_t_ps(const double t_ps) {
    t_ps_ = t_ps;                                                                  
  }

  double get_t_pp() const {
    return t_pp_;
  }
  void set_t_pp(const double t_pp) {
    t_pp_ = t_pp;
  }

  double get_t_ss100() const {
    return t_ss100_;
  }
  void set_t_ss100(const double t_ss100) {
    t_ss100_ = t_ss100;
  }

  double get_t_ss110() const {
    return t_ss110_;
  }
  void set_t_ss110(const double t_ss110) {
    t_ss110_ = t_ss110;
  }

  double get_t_ss101() const {
    return t_ss101_;
  }
  void set_t_ss101(const double t_ss101) {
    t_ss101_ = t_ss101;
  }

  double get_t_ss001() const {
    return t_ss001_;
  }
  void set_t_ss001(const double t_ss001) {
    t_ss001_ = t_ss001;
  }

  double get_t_ss111() const {
    return t_ss111_;
  }
  void set_t_ss111(const double t_ss111) {
    t_ss111_ = t_ss111;
  }

  double get_ep_d() const {
    return ep_d_;
  }
  void set_ep_d(const double ep_d) {
    ep_d_ = ep_d;
  }

  double get_es() const {
    return es_;
  }
  void set_es(const double es) {
    es_ = es;
  }

  double get_ep_p() const {
    return ep_p_;
  }
  void set_ep_p(const double ep_p) {
    ep_p_ = ep_p;
  }

  double get_U_dd() const {
    return U_dd_;
  }
  void set_U_dd(const double U_dd) {
    U_dd_ = U_dd;
  }

  double get_U_pp() const {
    return U_pp_;
  }
  void set_U_pp(const double U_pp) {
    U_pp_ = U_pp;
  }

private:
  double t_pd_ = 0;
  double t_ps_ = 0;
  double t_pp_ = 0;
  double t_ss100_ = 0;
  double t_ss110_ = 0;
  double t_ss101_ = 0;
  double t_ss001_ = 0;
  double t_ss111_ = 0;
  double ep_p_ = 0;
  double ep_d_ = 0;
  double es_ = 0;
  double U_dd_ = 0;
  double U_pp_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::dps_3D_model<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_pd_);
  buffer_size += concurrency.get_buffer_size(t_ps_);
  buffer_size += concurrency.get_buffer_size(t_pp_);
  buffer_size += concurrency.get_buffer_size(t_ss100_);
  buffer_size += concurrency.get_buffer_size(t_ss110_);
  buffer_size += concurrency.get_buffer_size(t_ss101_);
  buffer_size += concurrency.get_buffer_size(t_ss001_);
  buffer_size += concurrency.get_buffer_size(t_ss111_);
  buffer_size += concurrency.get_buffer_size(ep_p_);
  buffer_size += concurrency.get_buffer_size(ep_d_);
  buffer_size += concurrency.get_buffer_size(es_);
  buffer_size += concurrency.get_buffer_size(U_dd_);
  buffer_size += concurrency.get_buffer_size(U_pp_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::dps_3D_model<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_pd_);
  concurrency.pack(buffer, buffer_size, position, t_ps_);
  concurrency.pack(buffer, buffer_size, position, t_pp_);
  concurrency.pack(buffer, buffer_size, position, t_ss100_);
  concurrency.pack(buffer, buffer_size, position, t_ss110_);
  concurrency.pack(buffer, buffer_size, position, t_ss101_);
  concurrency.pack(buffer, buffer_size, position, t_ss001_);
  concurrency.pack(buffer, buffer_size, position, t_ss111_);
  concurrency.pack(buffer, buffer_size, position, ep_p_);
  concurrency.pack(buffer, buffer_size, position, ep_d_);
  concurrency.pack(buffer, buffer_size, position, es_);
  concurrency.pack(buffer, buffer_size, position, U_dd_);
  concurrency.pack(buffer, buffer_size, position, U_pp_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::dps_3D_model<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_pd_);
  concurrency.unpack(buffer, buffer_size, position, t_ps_);
  concurrency.unpack(buffer, buffer_size, position, t_pp_);
  concurrency.unpack(buffer, buffer_size, position, t_ss100_);
  concurrency.unpack(buffer, buffer_size, position, t_ss110_);
  concurrency.unpack(buffer, buffer_size, position, t_ss101_);
  concurrency.unpack(buffer, buffer_size, position, t_ss001_);
  concurrency.unpack(buffer, buffer_size, position, t_ss111_);
  concurrency.unpack(buffer, buffer_size, position, ep_p_);
  concurrency.unpack(buffer, buffer_size, position, ep_d_);
  concurrency.unpack(buffer, buffer_size, position, es_);
  concurrency.unpack(buffer, buffer_size, position, U_dd_);
  concurrency.unpack(buffer, buffer_size, position, U_pp_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::dps_3D_model<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("dps-3D-model");

    try {
      reader_or_writer.execute("t_pd", t_pd_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_ps", t_ps_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_pp", t_pp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_ss100", t_ss100_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_ss110", t_ss110_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_ss101", t_ss101_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_ss001", t_ss001_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_ss111", t_ss111_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_d", ep_d_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("es", es_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_p", ep_p_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U_dd", U_dd_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U_pp", U_pp_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
