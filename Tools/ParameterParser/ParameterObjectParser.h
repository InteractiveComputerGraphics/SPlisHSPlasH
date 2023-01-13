#pragma once

#include "SPlisHSPlasH/Common.h"
#include "ParameterObject.h"


namespace SPH
{
	class ParameterObjectParser
	{
	protected:
		typedef std::function<void(GenParam::NumericParameter<Real>*)>				RealParamCB;
		typedef std::function<void(GenParam::NumericParameter<unsigned int>*)>		UInt32ParamCB;
		typedef std::function<void(GenParam::NumericParameter<unsigned short>*)>	UInt16ParamCB;
		typedef std::function<void(GenParam::NumericParameter<unsigned char>*)>		UInt8ParamCB;
		typedef std::function<void(GenParam::NumericParameter<int>*)>				Int32ParamCB;
		typedef std::function<void(GenParam::NumericParameter<short>*)>				Int16ParamCB;
		typedef std::function<void(GenParam::NumericParameter<char>*)>				Int8ParamCB;
		typedef std::function<void(GenParam::EnumParameter*)>						EnumParamCB;
		typedef std::function<void(GenParam::BoolParameter*)>						BoolParamCB;
		typedef std::function<void(GenParam::StringParameter*)>						StringParamCB;
		typedef std::function<void(GenParam::RealVectorParameter*)>					VecRealParamCB;
		typedef std::function<void(GenParam::VectorParameter<unsigned int>*)>		VecUintParamCB;
		typedef std::function<void(const GenParam::	ParameterObject*)>				ParameterObjectCB;

		std::vector<RealParamCB> m_realParamCB;
		std::vector<UInt32ParamCB> m_uint32ParamCB;
		std::vector<UInt16ParamCB> m_uint16ParamCB;
		std::vector<UInt8ParamCB> m_uint8ParamCB;
		std::vector<Int32ParamCB> m_int32ParamCB;
		std::vector<Int16ParamCB> m_int16ParamCB;
		std::vector<Int8ParamCB> m_int8ParamCB;
		std::vector<EnumParamCB> m_enumParamCB;
		std::vector<BoolParamCB> m_boolParamCB;
		std::vector<StringParamCB> m_stringParamCB;
		std::vector<VecRealParamCB> m_vecRealParamCB;
		std::vector<VecUintParamCB> m_vecUintParamCB;
		std::vector<ParameterObjectCB> m_paramObjCB;

	public:
		void parseParameterObject(GenParam::ParameterObject* paramObj);
		void parseParameters();

		void addRealParamCB(RealParamCB cb) { m_realParamCB.push_back(cb); }
		void addUInt32ParamCB(UInt32ParamCB cb) { m_uint32ParamCB.push_back(cb); }
		void addUInt16ParamCB(UInt16ParamCB cb) { m_uint16ParamCB.push_back(cb); }
		void addUInt8ParamCB(UInt8ParamCB cb) { m_uint8ParamCB.push_back(cb); }
		void addInt32ParamCB(Int32ParamCB cb) { m_int32ParamCB.push_back(cb); }
		void addInt16ParamCB(Int16ParamCB cb) { m_int16ParamCB.push_back(cb); }
		void addInt8ParamCB(Int8ParamCB cb) { m_int8ParamCB.push_back(cb); }
		void addEnumParamCB(EnumParamCB cb) { m_enumParamCB.push_back(cb); }
		void addBoolParamCB(BoolParamCB cb) { m_boolParamCB.push_back(cb); }
		void addStringParamCB(StringParamCB cb) { m_stringParamCB.push_back(cb); }
		void addVecRealParamCB(VecRealParamCB cb) { m_vecRealParamCB.push_back(cb); }
		void addVecUintParamCB(VecUintParamCB cb) { m_vecUintParamCB.push_back(cb); }

		void addParameterObjectCB(ParameterObjectCB cb) { m_paramObjCB.push_back(cb); }
	};
}
