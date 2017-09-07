
#ifndef APICALFORCE_HPP_
#define APICALFORCE_HPP_

#include "AbstractForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law between two capsules (cylinder with hemispherical caps), defined in Farrell et al:
 * J R Soc Interface. 2017 Jun;14(131). pii: 20170073. doi: 10.1098/rsif.2017.0073
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class ApicalForce : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
    //friend class TestApicalForce;

private:

    /** The elastic modulus of both cells (Farrell et al) */
    double mSpringConstant;

    /** The elastic modulus of both cells (Farrell et al) */
    double mSpringLength;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mSpringConstant;
        archive & mSpringLength;

    }





public:

    /**
     * Constructor.
     */
    ApicalForce();

    /**
     * Destructor.
     */
    virtual ~ApicalForce() = default;

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

    void SetSpringConstant(double springConstant);
	double GetSpringConstant();
    void SetSpringLength(double springLength);
    double GetSpringLength();



};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ApicalForce)

#endif /*APICALFORCE_HPP_*/
