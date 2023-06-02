// Author: Clara Bayley
// File: "singlevarstorage.hpp"
/* Classes and helper functions in a namespace
useful for using storage clases with buffers to
write values of 1D data into chunks of arrays
in a zarr store */

#ifndef SINGLEVARSTORAGE_HPP
#define SINGLEVARSTORAGE_HPP

#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>
#include <algorithm>

#include "./zarrstores.hpp"

template <typename T>
class SingleVarStorage
{
private:
  virtual void writechunk() = 0;
  virtual void writejsons() = 0;

protected:
  FSStore &store;            // file system store satisfying zarr store specificaiton v2
  const std::string name;    // name to call variable being stored
  const std::string units;   // units of coordinate being stored (for arrayattrs json)
  const double scale_factor; // scale_factor of data (for array .zattrs json)
  std::vector<T> buffer;     // buffer to store values in until writing to array chunk

  const size_t chunksize;  // fixed size of array chunks (=max no. datapoints in buffer before writing)
  unsigned int chunkcount; // number of chunks of array so far written to store
  unsigned int bufferfill; // number of datapoints so far copied into buffer
  unsigned int ndata;      // number of data points that have been observed

  char zarr_format = '2';    // storage spec. version 2
  char order = 'C';                // layout of bytes within each chunk of array in storage, can be 'C' or 'F'
  std::string compressor = "null"; // compression of data when writing to store
  std::string fill_value = "null"; // fill value for empty datapoints in array
  std::string filters = "null";    // codec configurations for compression
  std::string dtype;               // datatype stored in arrays

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
  SingleVarStorage(FSStore &store, const unsigned int maxchunk,
                   const std::string name, const std::string dtype,
                   const std::string units, const double scale_factor)
      : store(store), name(name), units(units),
        scale_factor(scale_factor),
        buffer(maxchunk, std::numeric_limits<T>::max()),
        chunksize(maxchunk), chunkcount(0),
        bufferfill(0), ndata(0), dtype(dtype) {}

  virtual ~SingleVarStorage(){};

  void is_name(const std::string goodname) const
  {
    if (name != goodname)
    {
      const std::string errmsg("name of storage is " + name +
                               ", but should be " + goodname);
      throw std::invalid_argument(errmsg);
    }
  }

  int get_ndata() const { return ndata; }

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
    this->chunkcount = storagehelper::
        writebuffer2chunk(this->store, this->buffer,
                          this->name, this->chunkcount);
    this->bufferfill = 0;

    writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST) */
  {
    const auto shape("[" + std::to_string(this->ndata) + "]");
    const auto chunks("[" + std::to_string(this->chunksize) + "]");
    const std::string dims = "[\"" + this->name + "\"]";

    this->zarrayjsons(shape, chunks, dims);
  }

public:
  CoordinateStorage(FSStore &store, const unsigned int maxchunk,
                    const std::string name, const std::string dtype,
                    const std::string units, const double scale_factor)
      : SingleVarStorage<T>(store, maxchunk, name, dtype,
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
{
private:
  const unsigned int ngridboxes; // number of gridboxes that are observed (at each obs)

  void writechunk()
  /* write data in buffer to a chunk in store alongside metadata jsons */
  {
    const std::string chunknum = std::to_string(this->chunkcount) + ".0";
    this->chunkcount = storagehelper::
        writebuffer2chunk(this->store, this->buffer,
                          this->name, chunknum,
                          this->chunkcount);
    this->bufferfill = 0;

    writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST).
  Assert also check 2D data dimensions is as expected */
  {
    assert((this->ndata == nobs * ngridboxes) &&
           "1D data length matches 2D array size");
    assert((this->chunksize % ngridboxes == 0.0) &&
           "chunks are integer multple of number of gridboxes");

    const auto ngstr = std::to_string(ngridboxes);
    const auto nobstr = std::to_string(nobs);
    const auto nchstr = std::to_string(this->chunksize / ngridboxes);

    const auto shape("[" + nobstr + ", " + ngstr + "]");
    const auto chunks("[" + nchstr + ", " + ngstr + "]");
    const std::string dims = "[\"time\", \"gbxindex\"]";
    this->zarrayjsons(shape, chunks, dims);
  }

public:
  unsigned int nobs; // number of output times that have been observed

  TwoDStorage(FSStore &store, const unsigned int maxchunk,
              const std::string name, const std::string dtype,
              const std::string units, const double scale_factor,
              const unsigned int ngrid)
      : SingleVarStorage<T>(store, floor(maxchunk / ngrid) * ngrid,
                            name, dtype, units, scale_factor),
        ngridboxes(ngrid), nobs(0) {}

  ~TwoDStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      writechunk();
    }
  }
};

#endif // SINGLEVARSTORAGE_HPP