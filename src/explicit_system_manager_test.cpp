//
// Created by ryan on 3/9/16.
//

#include <gtest/gtest.h>
#include "explicit_system_manager.h"

using namespace explicit_fea;

TEST(ExplicitSystemManager, RunsSimulation) {
    std::string filename("/home/ryan/Dropbox/C++/explicit-beam-fea-v2/examples/single-element/config.json");
    ExplicitSystemManager esm(filename);
    esm.run();

	// Check displacements
	ColumnVector actual_displacements = esm.getExplicitSystem().getDisplacements();
	ColumnVector expected_displacements(ColumnVector::Zero(actual_displacements.size()));
	const std::vector<std::unique_ptr<BC>> &bcs = esm.getExplicitSystem().getMesh().get_bcs();
	expected_displacements[6] = esm.getIterationNumber() * esm.getTimeStep() * bcs[bcs.size() - 1]->get_value(esm.getExplicitSystem().getTime());
	for (int i = 0; i < expected_displacements.size(); ++i)
		EXPECT_FLOAT_EQ(expected_displacements[i], actual_displacements[i]);
}


