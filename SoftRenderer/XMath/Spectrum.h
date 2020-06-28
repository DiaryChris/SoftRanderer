#pragma once

#include <map>
#include "../XMath/Common.h"
#include "../XMath/XYZColorSpace.h"
#include "../XMath/ColorConversions.h"

namespace xmath
{
	using namespace std;

	class Spectrum
	{
	public:

		virtual ~Spectrum(){}

		virtual void addSample(float wavelen, float power) = 0;

		virtual float evaluate(float wavelen) const = 0;

		virtual void mul(Spectrum& r) = 0;

		// TODO IMPLEMENT!
		virtual float3 toXYZ()
		{
			float3 XYZ{0};

			static const float start = 360;
			static const float end = 830;
			static const float increment = 1;

			for (float wavelen = start; wavelen < end; wavelen += increment)
			{
				auto power = evaluate(wavelen);

				XYZ.x += xyz1931::X(wavelen) * power * increment;
				XYZ.y += xyz1931::Y(wavelen) * power * increment;
				XYZ.z += xyz1931::Z(wavelen) * power * increment;
			}

			return XYZ;
		}
	};

	// TODO 这里面的一些设计可以抽象出来，作为日后的 Track 的基础
	class SampledSpectrum:
		public Spectrum
	{
	public:
		SampledSpectrum(){}
		SampledSpectrum(map<float,float> samples):
			_samples(move(samples)){}
		SampledSpectrum(const SampledSpectrum& r):
			_samples(r._samples){}
		SampledSpectrum(SampledSpectrum&& r):
			_samples(move(r._samples)){}

		virtual void addSample(float wavelen, float power) override
		{
			_samples[wavelen] = power;
		}

		virtual float evaluate(float wavelen) const override
		{
			if (_samples.size() == 0)
				return 0;

			auto lowerBound = _samples.lower_bound(wavelen);

			if (lowerBound == _samples.end())
				return (--_samples.end())->second;

			auto prevSampleItr = lowerBound;
			prevSampleItr--;

			if (prevSampleItr == _samples.end())
				return lowerBound->second;

			auto t = (wavelen - prevSampleItr->first) / (lowerBound->first - prevSampleItr->first);

			return xmath::lerp(prevSampleItr->second, lowerBound->second, t);
		}

		virtual void mul(Spectrum& r) override
		{
			for (auto itr = _samples.begin(); itr != _samples.end(); ++itr)
			{
				itr->second *= r.evaluate(itr->first);
			}
		}

		map<float, float> _samples;
	};
}
