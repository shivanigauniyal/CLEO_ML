// Author: Clara Bayley
// File: "singlevarstorage.hpp"
/* Classes and helper functions
writing a single variable (could
be 1D or 2D, coordinate etc.)
via a buffer into chunks of arrays
in a zarr store */

#ifndef SINGLEVARSTORAGE_HPP
#define SINGLEVARSTORAGE_HPP

#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>
#include <algorithm>
#include <utility>
#include <tuple>
#include <stdexcept>

#include "./zarrstores.hpp"

inline unsigned int NOTSETCHUNKSIZE() { return std::numeric_limits<unsigned int>::min(); }
inline unsigned int NOTSETVALUE() { return std::numeric_limits<unsigned int>::max(); }

template <typename T>
class SingleVarStorage
{
private:
  virtual void writechunk() = 0;
  virtual void writejsons() = 0;

  size_t chunksize; // fixed size of array chunks (=max no. datapoints in buffer before writing)

protected:
  FSStore &store;            // file system store satisfying zarr store specificaiton v2
  const std::string name;    // name to call variable being stored
  const std::string units;   // units of coordinate being stored (for arrayattrs json)
  const double scale_factor; // scale_factor of data (for array .zattrs json)
  std::vector<T> buffer;     // buffer to store values in until writing to array chunk

  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed

  const char zarr_format = '2';          // storage spec. version 2
  const char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  const std::string compressor = "null"; // compression of data when writing to store
  const std::string fill_value = "null"; // fill value for empty datapoints in array
  const std::string filters = "null";    // codec configurations for compression
  const std::string dtype;               // datatype stored in arrays

  unsigned int get_chunksize() const { return chunksize; }

  void set_buffer_chunksize(const unsigned int i_chunksize)
  {
    if (chunksize != NOTSETCHUNKSIZE())
    {
      const std::string err("chunksize already set; it cannot be changed");
      throw std::invalid_argument(err);
    }
    chunksize = i_chunksize;
    buffer.assign(chunksize, std::numeric_limits<T>::max());
  }

  void zarrayjsons(const std::string shape,
                   const std::string chunks,
                   const std::string dims)
  /* write array's metadata to .json files */
  {
    const std::string metadata = storagehelper::
        metadata(zarr_format, order, shape, chunks, dtype,
                 compressor, fill_value, filters);

    const std::string arrayattrs = storagehelper::
        arrayattrs(dims, units, scale_factor);

    storagehelper::writezarrjsons(store, name, metadata, arrayattrs);
  }

  void copy2buffer(const T val)
  /* copy value 'val' to buffer */
  {
    bufferfill = storagehelper::val2buffer<T>(val, buffer, bufferfill);
    ++ndata;
  }

  void copy2buffer(const std::vector<T> &vec)
  /* copy values of type T in vector 'vec' to buffer */
  {
    bufferfill = storagehelper::vec2buffer<T>(vec, buffer, bufferfill);
    ndata += vec.size();
  }

public:
  SingleVarStorage(FSStore &store, const unsigned int chunksize,
                   const std::string name, const std::string dtype,
                   const std::string units, const double scale_factor)
      : chunksize(chunksize), store(store),
        name(name), units(units), scale_factor(scale_factor),
        buffer(chunksize, std::numeric_limits<T>::max()),
        chunkcount(0), bufferfill(0), ndata(0), dtype(dtype) {}

  virtual ~SingleVarStorage(){};

  int get_ndata() const { return ndata; }

  void is_name(const std::string &goodname) const
  {
    if (name != goodname)
    {
      const std::string errmsg("name of storage is " + name +
                               ", but should be " + goodname);
      throw std::invalid_argument(errmsg);
    }
  }

  void value_to_storage(const T val)
  /* write val in the zarr store. First copy it to a buffer,
  then write buffer to a chunk in the store when the number
  of values in the buffer reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      writechunk();
    }

    copy2buffer(val);
  }

  void value_to_storage(const std::vector<T> &vec)
  /* write 'vec' vector of type T in the zarr store.
  First copy vector to a buffer, then write buffer to a
  chunk in the store when the number of values in
  the buffer reaches the chunksize */
  {
    if (bufferfill == chunksize)
    {
      writechunk();
    }

    copy2buffer(vec);
  }
};

template <typename T>
struct CoordinateStorage : SingleVarStorage<T>
/* storage of a 1D variable with 'dims' in .zattrs metadata
equal to name of variable (ie. variable is an xarray coord)*/
{
private:
  void writechunk()
  /* write data in buffer to a chunk in store */
  {
    std::tie(this->chunkcount, this->bufferfill) =
        storagehelper::
            writebuffer2chunk(this->store, this->buffer,
                              this->name, this->chunkcount);

    writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST) */
  {
    const auto shape("[" + std::to_string(this->ndata) + "]");
    const auto chunks("[" + std::to_string(this->get_chunksize()) + "]");
    const std::string dims = "[\"" + this->name + "\"]";

    this->zarrayjsons(shape, chunks, dims);
  }

public:
  CoordinateStorage(FSStore &store, const unsigned int chunksize,
                    const std::string name, const std::string dtype,
                    const std::string units, const double scale_factor)
      : SingleVarStorage<T>(store, chunksize, name, dtype,
                            units, scale_factor) {}

  ~CoordinateStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      writechunk();
    }
  }
};

template <typename T>
struct TwoDStorage : SingleVarStorage<T>
/* 2D storage with dimensions [time, dim1] where
ntime is number of observation events (no. time outputs)
and ndim1 is the number of elements in 1st dimension
of 2-D data i.e. no. elements observed for each time.
For example, ndim1 could equal the number of gridboxes
an observer observes during 1 observation. Data for values
of time and dim1 could be output using a CoordinateStorage */
{
private:
  const std::string dim1name; // name of 1st dimension (e.g. "gbxindex")
  unsigned int ndim1;         // number elements in 1st dimensin (e.g. number of gridboxes that are observed)

  void writechunk()
  /* write data in buffer to a chunk in store alongside metadata jsons */
  {
    const std::string chunknum = std::to_string(this->chunkcount) + ".0";
    std::tie(this->chunkcount, this->bufferfill) = storagehelper::
        writebuffer2chunk(this->store, this->buffer,
                          this->name, chunknum,
                          this->chunkcount);

    writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST).
  Assert also check 2D data dimensions is as expected */
  {
    assert((this->ndata == nobs * ndim1) &&
           "1D data length matches 2D array size");
    assert((this->get_chunksize() % ndim1 == 0.0) &&
           "chunks are integer multiple of 1st dimension of 2-D data");

    const auto n1str = std::to_string(ndim1);
    const auto nobstr = std::to_string(nobs);
    const auto nchstr = std::to_string(this->get_chunksize() / ndim1);

    const auto shape("[" + nobstr + ", " + n1str + "]");
    const auto chunks("[" + nchstr + ", " + n1str + "]");
    const std::string dims = "[\"time\", \"" + dim1name + "\"]";
    this->zarrayjsons(shape, chunks, dims);
  }

protected:
  void set_ndim1(const unsigned int i_ndim1)
  {
    if (ndim1 != NOTSETVALUE())
    {
      const std::string err("ndim1 already set; it cannot be changed");
      throw std::invalid_argument(err);
    }
    ndim1 = i_ndim1;
  }

public:
  unsigned int nobs; // number of output times that have been observed

  TwoDStorage(FSStore &store, const unsigned int maxchunk,
              const std::string name, const std::string dtype,
              const std::string units, const double scale_factor,
              const std::string i_dim1name, const unsigned int i_ndim1)
      : SingleVarStorage<T>(store, storagehelper::good2Dchunk(maxchunk, i_ndim1),
                            name, dtype, units, scale_factor),
        dim1name(i_dim1name), ndim1(i_ndim1), nobs(0) {}

  ~TwoDStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      writechunk();
    }
  }

  void is_dim1(const size_t goodndim1,
               const std::string &goodname) const
  {
    if ((size_t)ndim1 != goodndim1)
    {
      const std::string errmsg("ndim1 is" +
                               std::to_string(ndim1) +
                               ", but should be " +
                               std::to_string(goodndim1));
      throw std::invalid_argument(errmsg);
    }

    if (dim1name != goodname)
    {
      const std::string errmsg("name of dim1 is " + dim1name +
                               ", but should be " + goodname);
      throw std::invalid_argument(errmsg);
    }
  }
};

CoordinateStorage<double> make_timezarr(FSStore &store,
                                        const int chunksize)
{
  const std::string name("time");
  const std::string dtype("<f8");
  const std::string units("s");
  constexpr double scale_factor(dlc::TIME0);

  return CoordinateStorage<double>(store, chunksize, name,
                                   dtype, units, scale_factor);
}

CoordinateStorage<unsigned int> make_gbxzarr(FSStore &store,
                                             const int chunksize)
{
  const std::string name("gbxindex");
  const std::string dtype("<u4");
  const std::string units(" ");
  constexpr double scale_factor(1);

  return CoordinateStorage<unsigned int>(store, chunksize, name,
                                         dtype, units, scale_factor);
}

TwoDStorage<size_t> make_nsuperszarr(FSStore &store,
                                     const int maxchunk,
                                     const unsigned int ngbxs)
{
  const std::string name("nsupers");
  const std::string dtype("<u8");
  const std::string units(" ");
  constexpr double scale_factor(1);
  const std::string dim1name("gbxindex");

  return TwoDStorage<size_t>(store, maxchunk, name, dtype, units,
                             scale_factor, dim1name, ngbxs);
}

TwoDStorage<size_t> make_nrainsuperszarr(FSStore &store,
                                     const int maxchunk,
                                     const unsigned int ngbxs)
{
  const std::string name("nrainsupers");
  const std::string dtype("<u8");
  const std::string units(" ");
  constexpr double scale_factor(1);
  const std::string dim1name("gbxindex");

  return TwoDStorage<size_t>(store, maxchunk, name, dtype, units,
                             scale_factor, dim1name, ngbxs);
}

#endif // SINGLEVARSTORAGE_HPP