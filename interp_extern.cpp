#undef NDEBUG // force assertions
#include <assert.h>
#include <iostream>

#include <Halide.h>

using namespace std;

template <typename T>
static inline int bsearch(const T *a, const int len, const T val) {
    int l = 0;
    int u = len;
    int m;
    while (l < u) {
        m = (l + u) / 2;
        if (a[m] < val) {
            l = m + 1;
        } else {
            u = m;
        }
    }
    return l;
}

int bsearch_lut_extern(halide_buffer_t *xs, halide_buffer_t *xp, int xp_extent,
                       halide_buffer_t *out) {
    assert(xs->dimensions == 2);
    assert(xp->dimensions == 1);
    assert(out->dimensions == 2);
    if (xs->host == nullptr || xp->host == nullptr) {
#if 0
        cout << "bsearch_lut_extern: bounds query" << endl;
#endif
        if (xs->is_bounds_query()) {
            // input and output shapes are the same
            for (int i = 0; i < xs->dimensions; i++) {
                xs->dim[i].min = out->dim[i].min;
                xs->dim[i].extent = out->dim[i].extent;
            }
        }
        if (xp->is_bounds_query()) {
            xp->dim[0].min = 0;
            xp->dim[0].extent = xp_extent;
        }
    } else {
        assert(xs->host);
        assert(xp->host);
        assert(out->host);
        assert(xs->type == halide_type_of<double>());
        assert(xp->type == halide_type_of<double>());
        assert(out->type == halide_type_of<int>());
#if 0
        cout << "bsearch_lut_extern: xs[x]: [" << xs->dim[0].min << ", "
             << xs->dim[0].min + xs->dim[0].extent << ")" << endl;
        cout << "bsearch_lut_extern: xs[y]: [" << xs->dim[1].min << ", "
             << xs->dim[1].min + xs->dim[1].extent << ")" << endl;
        cout << "bsearch_lut_extern: xp[x]: [" << xp->dim[0].min << ", "
             << xp->dim[0].min + xp->dim[0].extent << ")" << endl;
        cout << "bsearch_lut_extern: out[x]: [" << out->dim[0].min << ", "
             << out->dim[0].min + out->dim[0].extent << ")" << endl;
        cout << "bsearch_lut_extern: out[y]: [" << out->dim[1].min << ", "
             << out->dim[1].min + out->dim[1].extent << ")" << endl;
#endif
        const double *_xp = (const double *)xp->host;
        for (int y = out->dim[1].min; y < out->dim[1].min + out->dim[1].extent; y++) {
            for (int x = out->dim[0].min; x < out->dim[0].min + out->dim[0].extent; x++) {
                const int coords[2] = {x, y};
                const double _xs = *(const double *)xs->address_of(coords);
                *(int *)out->address_of(coords) = bsearch(_xp, xp->dim[0].extent, _xs);
            }
        }
    }
    return 0;
}

// xs: The x-coordinates at which to evaluate the interpolated values.
// xp: The x-coordinates of the data points, must be increasing.
// fp: The y-coordinates of the data points, same length as xp.
// This probably won't work if you tile across xp/fp extent - need entire length
int interp_extern(halide_buffer_t *xs, halide_buffer_t *xp, halide_buffer_t *fp,
                  int N_fft, halide_buffer_t *out) {
    assert(xs->dimensions == 2);
    assert(xp->dimensions == 1);
    assert(fp->dimensions == 2);
    assert(out->dimensions == 2);
    if (xs->host == nullptr || xp->host == nullptr || fp->host == nullptr) {
#if 0
        cout << "interp_extern: bounds query" << endl;
#endif
        if (xs->is_bounds_query()) {
            // input and output shapes are the same
            for (int i = 0; i < xs->dimensions; i++) {
                xs->dim[i].min = out->dim[i].min;
                xs->dim[i].extent = out->dim[i].extent;
            }
        }
        if (xp->is_bounds_query()) {
            xp->dim[0].min = 0;
            xp->dim[0].extent = N_fft;
        }
        if (fp->is_bounds_query()) {
            fp->dim[0].min = 0;
            fp->dim[0].extent = N_fft;
            fp->dim[1].min = 0;
            fp->dim[1].extent = out->dim[1].extent;
        }
    } else {
        assert(xs->host);
        assert(xp->host);
        assert(fp->host);
        assert(out->host);
        assert(xs->type == halide_type_of<double>());
        assert(xp->type == halide_type_of<double>());
        assert(fp->type == halide_type_of<double>());
        assert(out->type == halide_type_of<double>());
#if 0
        cout << "interp_extern: xs[x]: [" << xs->dim[0].min << ", "
             << xs->dim[0].min + xs->dim[0].extent << ")" << endl;
        cout << "interp_extern: xs[y]: [" << xs->dim[1].min << ", "
             << xs->dim[1].min + xs->dim[1].extent << ")" << endl;
        cout << "interp_extern: xp[x]: [" << xp->dim[0].min << ", "
             << xp->dim[0].min + xp->dim[0].extent << ")" << endl;
        cout << "interp_extern: fp[x]: [" << fp->dim[0].min << ", "
             << fp->dim[0].min + fp->dim[0].extent << ")" << endl;
        cout << "interp_extern: fp[y]: [" << fp->dim[1].min << ", "
             << fp->dim[1].min + fp->dim[1].extent << ")" << endl;
        cout << "interp_extern: out[x]: [" << out->dim[0].min << ", "
             << out->dim[0].min + out->dim[0].extent << ")" << endl;
        cout << "interp_extern: out[y]: [" << out->dim[1].min << ", "
             << out->dim[1].min + out->dim[1].extent << ")" << endl;
#endif
        for (int y = out->dim[1].min; y < out->dim[1].min + out->dim[1].extent; y++) {
            // must be able to search across entire xp
            const double *_xp = (const double *)xp->host;
            // must be able to index across entire row of fp
            const int coords_y[2] = {0, y};
            const double *_fp = (const double *)fp->address_of(coords_y);
            for (int x = out->dim[0].min; x < out->dim[0].min + out->dim[0].extent; x++) {
                int coords[2] = {x, y};
                const double _xs = *(double *)xs->address_of(coords);
                double *_out = (double *)out->address_of(coords);
                const int i = bsearch(_xp, xp->dim[0].extent, _xs);
                if (_xs <= _xp[0]) { //  implies i == 0
                    cout << "interp_extern: clamped low" << endl;
                    *_out = _fp[0];
                } else if (_xs >= _xp[xp->dim[0].extent - 1]) { // implies i == xp->dim[0].extent
                    cout << "interp_extern: clamped high" << endl;
                    *_out = _fp[fp->dim[0].extent - 1];
                } else {
                    // assert(_xs >= _xp[i - 1]);
                    // assert(_xp[i] > _xp[i - 1]);
                    // y = mx + b
                    *_out = ((_fp[i] - _fp[i - 1]) / (_xp[i] - _xp[i - 1])) * (_xs - _xp[i -1]) + _fp[i - 1];
                }
            }
        }
    }
    return 0;
}
