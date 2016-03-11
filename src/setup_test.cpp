
#include <gtest/gtest.h>
#include "csv_parser.h"
#include "setup.h"

using namespace explicit_fea;

namespace {
    void writeStringToTxt(std::string filename, std::string data) {
        std::ofstream output_file;
        output_file.open(filename);

        if (!output_file.is_open()) {
            std::cerr << "Error opening file" << filename << ".\n";
        }
        else {
            output_file << data;
            output_file.close();
        }
    }
}

TEST(Setup, CreatesConfigFromJSON) {
    std::string json = "{\"nodes\":\"nodes_file\"}\n";
    std::string filename = "CreatesCorrectConfig.json";
    writeStringToTxt(filename, json);

    rapidjson::Document doc = parseJSONConfig(filename);

    EXPECT_TRUE(doc.HasMember("nodes"));
    std::string nodes_file(doc["nodes"].GetString());
    EXPECT_EQ("nodes_file", nodes_file);

    if (std::remove(filename.c_str()) != 0) {
        std::cerr << "Error removing test csv file " << filename << ".\n";
    }
}

TEST(Setup, CreatesNodesFromJSON) {
    std::string nodes_file = "CreatesCorrectNodes.csv";
    std::string json("{\"nodes\":\"" + nodes_file + "\"}");
    rapidjson::Document doc;
    doc.Parse(json.c_str());
    std::vector<std::vector<double> > expected = {{1, 2, 3},
                                                  {4, 5, 6}};

    CSVParser csv;
    csv.write(nodes_file, expected, 1, ",");

    std::vector<Node> nodes = createNodeVecFromJSON(doc);

    for (size_t i = 0; i < nodes.size(); ++i) {
       for (size_t j = 0; j < nodes[i].size(); ++j) {
           EXPECT_EQ(expected[i][j], nodes[i][j]);
       }
    }
    if (std::remove(nodes_file.c_str()) != 0) {
       std::cerr << "Error removing test csv file " << nodes_file << ".\n";
    }
}

TEST(Setup, CreatesElemsFromJSON) {
    std::string elems_file = "CreatesElems_elems.csv";
    std::string props_file = "CreatesElems_props.csv";
    std::string json("{\"elems\":\"" + elems_file + "\",\"props\":\"" + props_file + "\"}");

    rapidjson::Document doc;
    doc.Parse(json.c_str());
    std::vector<std::vector<unsigned int> > expected_elems = {{1, 2},
                                                              {2, 3}};

    std::vector<std::vector<double> > expected_props = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                                                        {11, 12, 13, 14, 15, 16, 17, 18, 19, 20}};

    CSVParser csv;
    csv.write(elems_file, expected_elems, 0, ",");
    csv.write(props_file, expected_props, 1, ",");

    std::vector<std::unique_ptr<BeamElement>> elems = createElemVecFromJSON(doc);

    for (size_t i = 0; i < expected_elems.size(); ++i) {
        for (size_t j = 0; j < expected_elems[i].size(); ++j)
        {
            EXPECT_EQ(expected_elems[i][j], elems[i]->get_node_numbers()[j]);
        }
    }

    for (size_t i = 0; i < expected_props.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected_props[i][0], elems[i]->get_props().youngs_modulus);
        EXPECT_DOUBLE_EQ(expected_props[i][1], elems[i]->get_props().shear_modulus);
        EXPECT_DOUBLE_EQ(expected_props[i][2], elems[i]->get_props().area);
        EXPECT_DOUBLE_EQ(expected_props[i][3], elems[i]->get_props().Iz);
        EXPECT_DOUBLE_EQ(expected_props[i][4], elems[i]->get_props().Iy);
        EXPECT_DOUBLE_EQ(expected_props[i][5], elems[i]->get_props().J);
        EXPECT_DOUBLE_EQ(expected_props[i][6], elems[i]->get_props().density);
        EXPECT_DOUBLE_EQ(expected_props[i][7], elems[i]->get_props().normal_vec[0]);
        EXPECT_DOUBLE_EQ(expected_props[i][8], elems[i]->get_props().normal_vec[1]);
        EXPECT_DOUBLE_EQ(expected_props[i][9], elems[i]->get_props().normal_vec[2]);
    }
    if (std::remove(elems_file.c_str()) != 0) {
        std::cerr << "Error removing test csv file " << elems_file << ".\n";
    }
    if (std::remove(props_file.c_str()) != 0) {
        std::cerr << "Error removing test csv file " << props_file << ".\n";
    }
}

TEST(Setup, CreatesBCsFromJSON) {
    std::string bcs_file = "CreatesBCs.csv";
    std::string json("{\"bcs\":\"" + bcs_file + "\"}");
    std::vector<std::vector<double> > expected = {{10, 20, 30, 0},
                                                  {40, 50, 60, 1}};
    CSVParser csv;
    csv.write(bcs_file, expected, 1, ",");

    rapidjson::Document doc;
    doc.Parse(json.c_str());
    std::vector<std::unique_ptr<BC>> bcs = createBCVecFromJSON(doc);

    for (size_t i = 0; i < bcs.size(); ++i) {
        EXPECT_EQ((unsigned int) expected[i][0], bcs[i]->node);
        EXPECT_EQ((unsigned int) expected[i][1], bcs[i]->dof);
        EXPECT_DOUBLE_EQ(expected[i][2], bcs[i]->get_value(0.0));
        EXPECT_EQ((BC::Type)(unsigned int) expected[i][3], bcs[i]->type);
    }
    if (std::remove(bcs_file.c_str()) != 0) {
        std::cerr << "Error removing test csv file " << bcs_file << ".\n";
    }
}

TEST(Setup, CreatesForcesFromJSON) {
    std::string forces_file = "CreatesForces.csv";
    std::string json("{\"forces\":\"" + forces_file + "\"}");
    std::vector<std::vector<double> > expected = {{10, 20, 30},
                                                  {40, 50, 60}};
    CSVParser csv;
    csv.write(forces_file, expected, 1, ",");

    rapidjson::Document doc;
    doc.Parse(json.c_str());
    std::vector<std::unique_ptr<Force>> forces = createForceVecFromJSON(doc);

    for (size_t i = 0; i < forces.size(); ++i) {
        EXPECT_EQ((unsigned int) expected[i][0], forces[i]->node);
        EXPECT_EQ((unsigned int) expected[i][1], forces[i]->dof);
        EXPECT_DOUBLE_EQ(expected[i][2], forces[i]->get_value(0.0));
    }

    if (std::remove(forces_file.c_str()) != 0) {
        std::cerr << "Error removing test csv file " << forces_file << ".\n";
    }
}

TEST(Setup, CreatesColumnVectorFromJSON) {
    std::string col_vec_file = "CreatesColumnVector.csv";
    std::string json("{\"cvec\":\"" + col_vec_file + "\"}");

    rapidjson::Document doc;
    doc.Parse(json.c_str());

    std::vector<double> expected = {10, 20, 30, 40, 50, 60};
    std::string expected_string;
    for (size_t i = 0; i < expected.size(); ++i)
    {
		expected_string += std::to_string(expected[i]);
		if (i < expected.size() - 1)
			expected_string += ",";
    }
    writeStringToTxt(col_vec_file, expected_string);

    ColumnVector actual = createColumnVectorFromJSON(doc, "cvec", expected.size());

    for (size_t i = 0; i < expected.size(); ++i) {
       EXPECT_DOUBLE_EQ(expected[i], actual[i]);
    }
    if (std::remove(col_vec_file.c_str()) != 0) {
        std::cerr << "Error removing test csv file " << col_vec_file << ".\n";
    }
}