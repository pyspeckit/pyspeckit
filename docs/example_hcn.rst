Radio Fitting: HCN example with freely varying hyperfine amplitudes
===================================================================

.. include:: <isogrk3.txt>

Example hyperfine line fitting for the HCN 1-0 line.


.. literalinclude:: ../examples/hcn_example.py
   :language: python

The green lines in the following figures all show the residuals to the fit              

.. figure:: images/hcn_fixedhf_fit.png
        :alt: Fit to the 3 hyperfine components of HCN 1-0 simultaneously with fixed amplitudes.  The (0)'s indicate that this is the 0'th velocity component being fit (though that velocity corresponds to the 12-01 component of the line)
        :figwidth: 650
        :width: 650

        Fit to the 3 hyperfine components of HCN 1-0 simultaneously with fixed amplitudes.  The (0)'s indicate that this is the 0'th velocity component being fit (though that velocity corresponds to the 12-01 component of the line)

.. figure:: images/hcn_freehf_fit.png
        :alt: Fit to the 3 hyperfine components of HCN 1-0 simultaneously with freely varying amplitudes.  The (0)'s indicate that this is the 0'th velocity component being fit
        :figwidth: 650
        :width: 650

        Fit to the 3 hyperfine components of HCN 1-0 simultaneously with freely varying amplitudes.  The (0)'s indicate that this is the 0'th velocity component being fit

.. figure:: images/hcn_freehf_ampandwidth_fit.png
        :alt: Fit to the 3 hyperfine components of HCN 1-0 simultaneously.  The widths are allowed to vary in this example.
        :figwidth: 650
        :width: 650

        Fit to the 3 hyperfine components of HCN 1-0 simultaneously.  The widths are allowed to vary in this example.

.. figure:: images/hcn_fixedhf_fit_2components.png
        :alt: Two-component fit with amplitude & width fixed.
        :figwidth: 650
        :width: 650

        A two-component fit to the same spectrum.  It appears to be a much
        better fit, hinting that there are indeed two components (though the
        fit is probably not unique)
