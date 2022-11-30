#pragma once

#include <gsElasticity/gsIterative_multiStep.h>

namespace gismo
{
    template <class T>
    gsIterative_multiStep<T>::gsIterative_multiStep(gsBaseAssembler<T>& assembler_):gsIterative<T>(assembler_)
    {
        iStep = 1;
    }

    template <class T>
    void gsIterative_multiStep<T>::solve(index_t numSteps)
    {
        while (m_status == working && iStep <= numSteps)
        {
            if (!compute())
            {
                m_status = bad_solution;
                goto abort;
            }
            if (m_options.getInt("Verbosity") == solver_verbosity::all)
                gsInfo << status() << std::endl;
            if (residualNorm < m_options.getReal("AbsTol") ||
                updateNorm < m_options.getReal("AbsTol") ||
                residualNorm / initResidualNorm < m_options.getReal("RelTol") ||
                updateNorm / initUpdateNorm < m_options.getReal("RelTol"))
            {
                iStep++;
                if(iStep > numSteps)
                    m_status = converged;
                numIterations = 0;
                break;
            }
            else if (numIterations == m_options.getInt("MaxIters"))
                m_status = interrupted;
        }

    abort:;
        if (m_options.getInt("Verbosity") != solver_verbosity::none)
        {
            if (m_status != working)
                gsInfo << status() << std::endl;
            else
                gsInfo << "Step " << iStep - 1 << " converged!" << std::endl;
        }
    }
}