"""
brukerread.py
=============
Python port of the KAZAN dataviewer Bruker EPR file readers.

Originally by Boris Epel & Alexey Silakov,
MPI of Bioinorganic Chemistry, Muelheim an der Ruhr, 2003-2004.

Public functions
----------------
xeprdsc(filename)        -> dsc  (dict of str->str)
xeprpar(dsc)             -> ax   (dict)
xeprparjss(dsc)          -> ax   (dict, JSS-only variant)
brukerread(filename)     -> (ax, y, dsc)  general reader
brukerreadjss(filename)  -> (ax, y, dsc)  JSS-only reader

All functions accept 1-3 return values just like their MATLAB counterparts,
but in Python we always return the full triple (ax, y, dsc) and let the
caller unpack what it needs.

ax  : dict with keys  x, xlabel, y, ylabel, z, zlabel, step, complex,
                       freq1 (optional), cf (optional), title (optional)
y   : numpy array, 1-D or 2-D
dsc : dict  str -> str  (raw parameters from .dsc / .par file)

Converted from Kazan Viewer with Claude Sonnet 4.5
Alexey Silakov 2026
"""

import os
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safeget(dsc, key, default=''):
    """Return dsc[key] if present, otherwise default."""
    return dsc.get(key, default)


def _stripunit(s):
    """Keep only digit / dot / E characters (mirrors MATLAB stripunit)."""
    return ''.join(c for c in s if c.isdigit() or c in '.E')


def _trim(s):
    return s.strip()


# ---------------------------------------------------------------------------
# xeprdsc  –  read .dsc / .par parameter file into a dict
# ---------------------------------------------------------------------------

def xeprdsc(filename):
    """
    Read a Bruker .dsc or .par file and return a dict of {key: value} strings.

    Parameters
    ----------
    filename : str

    Returns
    -------
    dsc : dict  {str: str}
    """
    forbidden = set('~!@#$%^&*()/\\.')
    dsc = {}
    last_key = ''

    with open(filename, 'r', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            parts = line.split(None, 1)  # split on first whitespace
            if not parts:
                continue
            a = parts[0]
            rest = parts[1] if len(parts) > 1 else ''

            if not a:
                continue
            if a[0] in forbidden:
                continue

            try:
                if not a[0].isdigit():
                    last_key = a
                    s = _trim(rest)
                    # strip surrounding single quotes
                    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
                        s = s[1:-1]
                    dsc[a] = s
                else:
                    # continuation line: append to the last key's value
                    if last_key and last_key in dsc:
                        dsc[last_key] = _trim(dsc[last_key] + rest)
            except Exception:
                pass

    return dsc


# ---------------------------------------------------------------------------
# xeprpar  –  build axis dict from a dsc dict (general format)
# ---------------------------------------------------------------------------

def xeprpar(dsc):
    """
    Extract axis information from a dsc dict (general BES3T / ESP format).

    Parameters
    ----------
    dsc : dict

    Returns
    -------
    ax : dict
    """
    if not isinstance(dsc, dict):
        raise TypeError('dsc must be a dict')

    ax = {
        'x': np.array([]), 'xlabel': '',
        'y': np.array([]), 'ylabel': '',
        'z': np.array([]), 'zlabel': '',
        'step': 0,
        'complex': False,
    }

    # Frequency
    if 'MWFQ' in dsc:
        ax['freq1'] = float(dsc['MWFQ'])
    if 'MF' in dsc:
        ax['freq1'] = float(dsc['MF'])

    # Center field
    if 'CenterField' in dsc:
        ax['cf'] = float(_stripunit(dsc['CenterField'])) * 1e-4

    # Complex flag
    if 'IKKF' in dsc:
        ax['complex'] = (dsc['IKKF'] == 'CPLX')
    elif 'XQAD' in dsc:
        ax['complex'] = (dsc['XQAD'] == 'ON')

    # Title
    if 'TITL' in dsc:
        ax['title'] = dsc['TITL']
    elif 'JCO' in dsc:
        ax['title'] = dsc['JCO']
    else:
        ax['title'] = '?'

    ax = _getrespar321(dsc, ax, 'x')
    ax = _getrespar321(dsc, ax, 'y')

    return ax


def _getrespar321(dsc, ax, axletter):
    """Populate one axis (x or y) in ax dict from dsc parameters."""
    lo = axletter.lower()
    up = axletter.upper()

    ax[lo] = np.array([])
    ax[lo + 'label'] = ''

    dim = 1
    sf = 0.0
    step = 1.0
    label = '?'

    unit = _safeget(dsc, up + 'UNI', '?').replace("'", '')

    # XSophe compatibility shim
    if 'JSS' not in dsc and ('RES' in dsc or 'RRES' in dsc) and 'HCF' in dsc:
        dsc = dict(dsc)  # shallow copy so we don't mutate the caller's dict
        dsc['JSS'] = '2'
        dsc['RES'] = _safeget(dsc, 'RES', _safeget(dsc, 'RRES', '1024'))

    if 'JSS' in dsc:
        jss = int(float(dsc['JSS']))

        if jss == 2:
            # CW field sweep
            if up != 'Y':
                dim = int(_safeget(dsc, 'RES', '1024'))
                if 'GST' in dsc and 'GSI' in dsc:
                    sf = float(dsc['GST'])
                    step = float(dsc['GSI']) / (dim - 1)
                elif 'HCF' in dsc and 'HSW' in dsc:
                    center = float(dsc['HCF'])
                    width = float(dsc['HSW'])
                    step = width / (dim - 1)
                    sf = center - width / 2
                elif 'HCF' in dsc and 'GST' in dsc:
                    center = float(dsc['HCF'])
                    width = 2 * abs(center - float(dsc['GST']))
                    step = width / (dim - 1)
                    sf = center - width / 2
                elif 'HCF' in dsc and 'GSI' in dsc:
                    dsc['DOS'] = '1'
                    center = float(dsc['HCF'])
                    width = float(dsc['GSI'])
                    step = width / (dim - 1)
                    sf = center - width / 2
                else:
                    # fallback – centre not defined, treat sf=0
                    step = 1.0
                    sf = 0.0
                label = 'Magnetic Field, G'

        elif jss == 32:
            # ESP 380 pulse
            if up != 'Y':
                if up + 'QNT' in dsc:
                    dim = int(dsc[up + 'PLS'])
                    idx = dsc[up + 'QNT'].replace(' ', '')
                    if idx == 'Time':
                        step = 0
                        val = np.fromstring(dsc.get('Psd5', ''), sep=' ')
                        for i in range(68, 75):
                            if i < len(val) and val[i] > 0:
                                ax['step'] = val[i]
                                break
                    elif idx == 'Magn.Field':
                        if 'HCF' in dsc:
                            cf = float(dsc['HCF'])
                            wd = float(dsc['HSW']) if 'HSW' in dsc else abs(cf - float(dsc.get('GST', str(cf)))) * 2
                        else:
                            cf = float(dsc['GST'])
                            wd = float(dsc['GSI'])
                            cf = cf + wd / 2
                        step = wd / (dim - 1)
                        sf = cf - wd / 2
                        label = 'Magnetic Field, G'
                    elif idx in ('RF1', 'RF2'):
                        sf = float(dsc.get(idx + 'StartFreq', '0').split()[0])
                        wd = float(dsc.get(idx + 'SweepWidth', '0').split()[0])
                        step = wd / (dim - 1)
                        label = idx + ', MHz'
                    elif idx == '1.RFSource':
                        sf = float(dsc.get('ESF', '0'))
                        wd = float(dsc.get('ESW', '0'))
                        step = wd / (dim - 1)
                        label = 'RF, MHz'
                else:
                    dim = int(dsc[up + 'PLS'])
                    if 'JUN' in dsc:
                        unit = dsc['JUN']
                    if 'GST' in dsc:
                        sf = float(dsc['GST'])
                        wd = float(dsc['GSI'])
                        step = wd / (dim - 1)
                        label = '?,' + unit
                    else:
                        if 'HCF' in dsc:
                            ax['cf'] = float(dsc['HCF']) * 1e-4
                        val = np.fromstring(dsc.get('XPD9', ''), sep=' ')
                        step = val[5] * 8 if len(val) > 5 else 1.0
                        label = 'Time, ns'

        elif jss in (4128, 4144):
            # 2-D files
            dim = int(dsc.get('SS' + up, '1'))
            if up == 'X' and ax.get('complex'):
                dim = dim // 2
            sw = float(_safeget(dsc, 'X' + up + 'WI', str(dim)))
            step = sw / (dim - 1) if dim > 1 else sw
            unit = _safeget(dsc, 'X' + up + 'UN', '?').replace("'", '')
            label = 'Time, ' + unit

        else:
            # Latest Bruker format / predefined ESP580 experiments
            if up + 'TYP' in dsc:
                if dsc[up + 'TYP'] != 'NODATA':
                    nam = _safeget(dsc, up + 'NAM', '?')
                    dim = int(dsc[up + 'PTS'])
                    sf = float(dsc[up + 'MIN'])
                    wd = float(dsc[up + 'WID'])
                    step = wd / (dim - 1)
                    label = nam + ' ' + unit
            elif up + 'NAM' in dsc:
                dim = int(dsc[up + 'PTS'])
                sf = float(dsc[up + 'MIN'])
                wd = float(dsc[up + 'WID'])
                step = wd / (dim - 1)
                nam = _safeget(dsc, up + 'NAM', '?')
                label = nam + ', ' + unit

    elif 'JEX' in dsc:
        if up != 'Y':
            dim = int(_safeget(dsc, 'RES', '1024'))
            if dsc['JEX'] == 'ENDOR':
                sf = float(dsc['ESF'])
                width = float(dsc['ESW'])
                step = width / (dim - 1)
                label = 'Frequency, MHz'

    elif up + 'TYP' in dsc:
        if dsc[up + 'TYP'] != 'NODATA':
            dim = int(dsc[up + 'PTS'])
            sf = float(dsc[up + 'MIN'])
            wd = float(dsc[up + 'WID'])
            step = wd / (dim - 1)
            nam = _safeget(dsc, up + 'NAM', '?')
            label = nam + ', ' + unit
            if up + 'AxisQuant' in dsc:
                idx = dsc[up + 'AxisQuant'].replace(' ', '')
                label = idx + ', ' + unit

    ax[lo] = sf + step * np.arange(dim)
    ax[lo + 'label'] = label
    return ax


# ---------------------------------------------------------------------------
# xeprparjss  –  JSS-only axis extractor
# ---------------------------------------------------------------------------

def xeprparjss(dsc):
    """
    Extract axis information from a dsc dict (JSS parameter variant).

    Parameters
    ----------
    dsc : dict

    Returns
    -------
    ax : dict
    """
    ax = {
        'x': np.array([]), 'xlabel': '',
        'y': np.array([]), 'ylabel': '',
        'z': np.array([]), 'zlabel': '',
        'step': 0,
    }

    if 'JSS' not in dsc:
        return ax

    jss = int(float(dsc['JSS']))

    if jss == 2:
        # CW field sweep
        dim = int(dsc['RES']) if 'RES' in dsc else 1024
        center = float(dsc['HCF'])
        width = float(dsc['HSW'])
        step = width / (dim - 1)
        ax['step'] = step
        ax['x'] = np.linspace(center - width / 2, center + width / 2, dim)
        ax['xlabel'] = 'Magnetic field, Gs'

    elif jss == 32:
        # ESP 380 pulse
        dim = int(dsc['XPLS'])
        val = np.fromstring(dsc.get('XPD9', ''), sep=' ')
        step = val[5] * 8 if len(val) > 5 else 1.0
        ax['step'] = step
        ax['x'] = step * np.arange(dim)
        ax['xlabel'] = 'Time, ns'

    else:
        # Generic fallback
        if 'RES' in dsc:
            dim = int(dsc['RES'])
        elif 'FRES' in dsc:
            dim = int(dsc['FRES'])
        elif 'XPLS' in dsc:
            dim = int(dsc['XPLS'])
        else:
            dim = 1024

        if 'GST' in dsc:
            cf = float(dsc['GST'])
            wd = float(dsc['GSI'])
            cf = cf + wd / 2
        else:
            cf = float(dsc['HCF'])
            if 'HSW' in dsc:
                wd = float(dsc['HSW'])
            else:
                wd = abs(cf - float(dsc.get('GST', str(cf)))) * 2

        step = wd / (dim - 1)
        sf = cf - wd / 2
        ax['step'] = step
        ax['x'] = sf + step * np.arange(dim)

        xlabel = dsc.get('JEX', '?')
        if 'JUN' in dsc:
            xlabel = xlabel + ', ' + dsc['JUN']
        else:
            xlabel = xlabel + ', ?'
        ax['xlabel'] = xlabel

    if 'MF' in dsc:
        ax['freq1'] = float(dsc['MF']) * 1e9

    return ax


# ---------------------------------------------------------------------------
# _getmatrix  –  read binary data array from file
# ---------------------------------------------------------------------------

def _getmatrix(filename, dims, fmt, byte_order, is_complex):
    """
    Read a flat binary file and reshape it to dims.

    Parameters
    ----------
    filename   : str
    dims       : tuple/list of int  (n_x, n_y, 1)
    fmt        : numpy dtype string  e.g. '<f4', '>i4', '<f8'
    byte_order : str  'ieee-le' or 'ieee-be'
    is_complex : bool

    Returns
    -------
    y : numpy array, shape (dims[0], dims[1]) or (dims[0],) if dims[1]==1
    """
    n_elements = int(np.prod(dims)) * (2 if is_complex else 1)
    raw = np.fromfile(filename, dtype=fmt, count=n_elements)

    if len(raw) < n_elements:
        raise IOError('Unable to read all expected data from ' + filename)

    if is_complex:
        raw = raw[0::2] + 1j * raw[1::2]

    # reshape; squeeze trailing size-1 dimensions
    total = int(np.prod(dims))
    if len(raw) >= total:
        raw = raw[:total]
    out = raw.reshape(dims[0], dims[1], order='F')
    if out.shape[1] == 1:
        out = out[:, 0]
    return out


def _numpy_fmt(fmt_str, byte_order):
    """Convert MATLAB-style format + byte_order to a numpy dtype."""
    endian = '<' if byte_order == 'ieee-le' else '>'
    mapping = {
        'float64': endian + 'f8',
        'float':   endian + 'f4',
        'float32': endian + 'f4',
        'int32':   endian + 'i4',
        'long':    endian + 'i4',
    }
    return mapping.get(fmt_str, endian + 'f8')


# ---------------------------------------------------------------------------
# brukerread  –  general reader
# ---------------------------------------------------------------------------

def brukerread(filename):
    """
    Read Bruker EPR data from .dsc/.dta or .par/.spc files (general format).

    Parameters
    ----------
    filename : str  (any of the four extensions)

    Returns
    -------
    ax  : dict    axis / metadata dictionary
    y   : ndarray 1-D or 2-D data array
    dsc : dict    raw parameter strings
    """
    ppath, name, ext = _splitext(filename)
    ext_up = ext.upper()

    if ext_up in ('.DSC', '.DTA'):
        data_file = os.path.join(ppath, name + '.dta')
        dsc_file  = os.path.join(ppath, name + '.dsc')
        file_type = 'BES3T'
    elif ext_up in ('.PAR', '.SPC'):
        data_file = os.path.join(ppath, name + '.spc')
        dsc_file  = os.path.join(ppath, name + '.par')
        file_type = 'PAR'
    else:
        raise ValueError('Unrecognised extension: ' + ext)

    dsc = xeprdsc(dsc_file)
    ax  = xeprpar(dsc)

    dims = [max(ax['x'].size, 1), max(ax['y'].size, 1), 1]
    byte_order = 'ieee-le'

    if file_type == 'BES3T':
        if _safeget(dsc, 'BSEQ', 'BIG') == 'BIG':
            byte_order = 'ieee-be'
        irfmt = _safeget(dsc, 'IRFMT', 'D')
        fmt = 'int32' if irfmt == 'I' else 'float64'
        np_fmt = _numpy_fmt(fmt, byte_order)
        y = _getmatrix(data_file, dims, np_fmt, byte_order, ax['complex'])

    else:  # PAR/SPC
        jss = int(float(_safeget(dsc, 'JSS', '2')))
        if jss > 0:
            dos_val = _safeget(dsc, 'DOS', '0')
            # DOS == 'Format' is the JSS indicator; numeric non-zero means big-endian
            try:
                if int(dos_val):
                    byte_order = 'ieee-be'
            except ValueError:
                pass  # 'Format' or other strings leave byte_order as-is

            if ax['complex']:
                # workaround: complex stored as interleaved real blocks
                dims2 = list(dims)
                dims2[0] *= 2
                np_fmt = _numpy_fmt('int32', byte_order)
                y2 = _getmatrix(data_file, dims2, np_fmt, byte_order, False)
                sz = ax['x'].size
                y = y2[:sz] + 1j * y2[sz:]
            else:
                vers = int(float(_safeget(dsc, 'VERS', '0')))
                fmt = 'float' if vers == 769 else 'float'
                np_fmt = _numpy_fmt(fmt, byte_order)
                y = _getmatrix(data_file, dims, np_fmt, byte_order, False)
        else:
            y = np.array([])

    # Non-uniformly spaced axes (BES3T IGD format)
    for axletter in ('X', 'Y'):
        lo = axletter.lower()
        if dsc.get(axletter + 'TYP') == 'IGD':
            try:
                igd_file = os.path.join(ppath, name + '.' + axletter + 'GF')
                tmp = np.fromfile(igd_file, dtype=(('>' if byte_order == 'ieee-be' else '<') + 'f8'))
                ax[lo] = tmp
            except Exception:
                print('Error reading axis file for', axletter)

    # Safety fallbacks
    if ax['x'].size == 0 or ax['x'].size != (y.shape[0] if y.ndim >= 1 else 0):
        ax['x'] = np.arange(1, (y.shape[0] if y.ndim >= 1 else 1) + 1, dtype=float)
    if ax['y'].size == 0 or ax['y'].size != (y.shape[1] if y.ndim == 2 else 0):
        ax['y'] = np.arange(1, (y.shape[1] if y.ndim == 2 else 1) + 1, dtype=float)

    return ax, y, dsc


# ---------------------------------------------------------------------------
# brukerreadjss  –  JSS-only reader
# ---------------------------------------------------------------------------

def brukerreadjss(filename):
    """
    Read Bruker EPR data from .dsc/.dta or .par/.spc files (JSS variant).

    Parameters
    ----------
    filename : str

    Returns
    -------
    ax  : dict    axis / metadata dictionary
    y   : ndarray 1-D or 2-D data array
    dsc : dict    raw parameter strings
    """
    ppath, name, ext = _splitext(filename)
    ext_up = ext.upper()

    if ext_up in ('.DSC', '.DTA'):
        data_file = os.path.join(ppath, name + '.dta')
        dsc_file  = os.path.join(ppath, name + '.dsc')
    elif ext_up in ('.PAR', '.SPC'):
        data_file = os.path.join(ppath, name + '.spc')
        dsc_file  = os.path.join(ppath, name + '.par')
    else:
        raise ValueError('Unrecognised extension: ' + ext)

    dsc = xeprdsc(dsc_file)
    ax  = xeprparjss(dsc)

    if 'DOS' not in dsc:
        dsc['DOS'] = '0'

    dos = dsc['DOS']
    if dos == 'Format':
        # Little-endian float32
        y = np.fromfile(data_file, dtype='<f4')
    else:
        # Big-endian int32 (long)
        y = np.fromfile(data_file, dtype='>i4')

    # Reshape to 2-D if y-axis has more than one point
    n_x = ax['x'].size if ax['x'].size > 0 else y.size
    n_y = ax['y'].size if ax['y'].size > 0 else 1

    if n_y > 1 and y.size == n_x * n_y:
        y = y.reshape(n_x, n_y, order='F')

    # Safety fallbacks
    n_rows = y.shape[0]
    n_cols = y.shape[1] if y.ndim == 2 else 1

    if ax['x'].size == 0 or ax['x'].size != n_rows:
        ax['x'] = np.arange(1, n_rows + 1, dtype=float)
    if ax['y'].size == 0 or ax['y'].size != n_cols:
        ax['y'] = np.arange(1, n_cols + 1, dtype=float)

    return ax, y, dsc


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _splitext(filename):
    """Return (directory, stem, extension) — mirrors MATLAB fileparts."""
    ppath = os.path.dirname(filename)
    base  = os.path.basename(filename)
    name, ext = os.path.splitext(base)
    return ppath, name, ext


# ---------------------------------------------------------------------------
# Quick self-test / usage example
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: python brukerread.py <filename.[dsc|dta|par|spc]>')
        sys.exit(0)

    fname = sys.argv[1]
    ax, y, dsc = brukerread(fname)
    print('Shape of y :', y.shape)
    print('x axis     :', ax['x'][:5], '...')
    print('x label    :', ax.get('xlabel', ''))
    print('Parameters :', list(dsc.keys())[:10], '...')