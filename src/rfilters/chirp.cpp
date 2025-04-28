#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _rfilter-chirp:

Chirp filter (:monosp:`chirp`)
----------------------------------

.. pluginparameters::

 * - lobes
   - |int|
   - Sets the desired number of filter side-lobes. The higher, the closer the
     filter will approximate an optimal low-pass filter, but this also increases
     ringing. Values of 2 or 3 are common (Default: 3)

This is a windowed version of the theoretically optimal low-pass filter. It
is generally one of the best available filters in terms of producing sharp
high-quality output. Its main disadvantage is that it produces strong ringing
around discontinuities, which can become a serious problem when rendering
bright objects with sharp edges (a directly visible light source will for
instance have black fringing artifacts around it). This is also the
computationally slowest reconstruction filter.

.. tabs::
    .. code-tab:: xml
        :name: lanczos-rfilter

        <rfilter type="lanczos">
            <integer name="lobes" value="4"/>
        </rfilter>

    .. code-tab:: python

        'type': 'lanczos',
        'lobes': 4

 */

template <typename Float, typename Spectrum>
class ChirpFilter final : public ReconstructionFilter<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ReconstructionFilter, init_discretization, m_radius)
    MI_IMPORT_TYPES()

    ChirpFilter(const Properties &props) : Base(props) {
        m_num_cycles = (ScalarFloat) props.get<ScalarFloat>("num_cycles", 3.0);
        m_central_frequency = (ScalarFloat) props.get<ScalarFloat>("central_frequency", 5000000.0);
        m_sampling_frequency = (ScalarFloat) props.get<ScalarFloat>("sampling_frequency", 10000000.0);
        m_speed_of_sound = (ScalarFloat) props.get<ScalarFloat>("speed_of_sound", 1540.0);

        m_radius = dr::ceil((m_sampling_frequency / m_central_frequency) * (m_num_cycles/2));
        // m_radius = m_num_cycles;

        init_discretization();
    }

    Float eval(Float x, dr::mask_t<Float> /*active*/ ) const override {
        
        Float x_signal = (m_radius + x) / m_sampling_frequency;
        Float x_gaussian = (x / m_radius) * 3;
        Float gaussian_part = dr::exp(-dr::sqr(x_gaussian));
        Float signal_part = dr::sin(2 * dr::Pi<Float> * m_central_frequency * x_signal);

        Float result = gaussian_part * signal_part;

        // self.ax_signal = lambda y: torch.exp(-1/2 * (3*(1/(self.tone_length/2))*y) ** 2) * torch.sin(2 * torch.pi * self.signal_frequency * y /  self.sound_speed)

        return dr::select(dr::abs(x) > m_radius, 0.f, result);
    }

    std::string to_string() const override {
        return tfm::format("ChirpFilter[m_radius=%f]", m_radius);
    }

    MI_DECLARE_CLASS()

    protected:
    ScalarFloat m_num_cycles;
    ScalarFloat m_central_frequency;
    ScalarFloat m_sampling_frequency;
    ScalarFloat m_speed_of_sound;
    Float m_tone_length;

};

MI_IMPLEMENT_CLASS_VARIANT(ChirpFilter, ReconstructionFilter)
MI_EXPORT_PLUGIN(ChirpFilter, "Chirp filter")
NAMESPACE_END(mitsuba)
