#pragma once

#if USE_DOUBLE
	#define CUNSEARCH_USE_DOUBLE_PRECISION
#endif

#if defined(USE_cuNSearch)
	#include "cuNSearch.h"
	typedef cuNSearch::NeighborhoodSearch NeighborhoodSearch;
#elif defined(USE_CompactNSearch)
	#include "CompactNSearch.h"
	typedef CompactNSearch::NeighborhoodSearch NeighborhoodSearch;
#elif defined(USE_TreeNSearch)
	#include "TreeNSearch.h"
	typedef tns::TreeNSearch NeighborhoodSearch;
#endif


namespace SPH
{
	/** This class is a wrapper for different neighborhood search methods.
	* The wrapper defines the interface of the neighborhood search.
	*/
	class NeighborhoodSearchWrapper
	{
	protected: 
		NeighborhoodSearch* m_neighborhoodSearch;
		std::vector<void*> m_mapPointSet_UserData;

	public:
		/** Creates and initializes a new instance of NeighborhoodSearch.
		* The caller has to ensure that this instance is deleted.
		*/
		NeighborhoodSearchWrapper(const Real supportRadius)
		{
			m_neighborhoodSearch = nullptr;

#if defined(USE_cuNSearch)
			m_neighborhoodSearch = new NeighborhoodSearch(supportRadius);
			m_neighborhoodSearch->set_radius(supportRadius);
#elif defined(USE_CompactNSearch)
			m_neighborhoodSearch = new NeighborhoodSearch(supportRadius, false);
			m_neighborhoodSearch->set_radius(supportRadius);
#elif defined(USE_TreeNSearch)
			m_neighborhoodSearch = new NeighborhoodSearch();
			const int n_threads = omp_get_num_procs();
			m_neighborhoodSearch->set_n_threads(n_threads);
			m_neighborhoodSearch->set_cell_size(1.5f * supportRadius);
			m_neighborhoodSearch->set_recursion_cap(1000);
			m_neighborhoodSearch->set_search_radius(supportRadius);
#endif
		}

		~NeighborhoodSearchWrapper()
		{
			delete m_neighborhoodSearch;
		}

		const std::vector<void*>& getMapPointSet_UserData() const
		{
			return m_mapPointSet_UserData;
		}

		void setNeighborhoodSearchActive(const bool active)
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->set_all_searches(active);
#else		
			m_neighborhoodSearch->set_active(active);
#endif
		}

		void setNeighborhoodSearchActive(unsigned int i, unsigned int j, bool active)
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->set_active_search((int)i, (int)j, active);
#else		
			m_neighborhoodSearch->set_active(i, j, active);
#endif
		}

		void setNeighborhoodSearchActive(unsigned int i, bool search_neighbors = true, bool find_neighbors = true)
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->set_active_search((int)i, search_neighbors, find_neighbors);
#else		
			m_neighborhoodSearch->set_active(i, search_neighbors, find_neighbors);
#endif
		}

		void zSort()
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->prepare_zsort();
#else		
			m_neighborhoodSearch->z_sort();
#endif
		}

		template <typename T>
		void applyZSort(const unsigned int pointSetIndex, T* data_ptr)
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->apply_zsort(pointSetIndex, data_ptr);
#else		
			auto const& d = m_neighborhoodSearch->point_set(pointSetIndex);
			d.sort_field(data_ptr);
#endif
		}

		void updatePointSets()
		{
#if defined(USE_TreeNSearch)
#else
			m_neighborhoodSearch->update_point_sets();
#endif
		}

		void findNeighbors()
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->run();
#else		
			m_neighborhoodSearch->find_neighbors();
#endif
		}

#if defined(USE_CompactNSearch)
		void find_neighbors(unsigned int point_set_id, unsigned int point_index, std::vector<std::vector<unsigned int>>& neighbors)
		{
			m_neighborhoodSearch->find_neighbors(point_set_id, point_index, neighbors);
		}

		void find_neighbors(Real const* x, std::vector<std::vector<unsigned int>>& neighbors)
		{
			m_neighborhoodSearch->find_neighbors(x, neighbors);
		}
#endif 

		unsigned int addPointSet(Real* x, std::size_t n, bool is_dynamic, bool search_neighbors, bool find_neighbors, void* user_data = nullptr)
		{
			unsigned int id = 0;
#if defined(USE_TreeNSearch)
			//if (is_dynamic) {
				id = m_neighborhoodSearch->add_point_set(x, (int)n);
			//}
			//else 
			//{
			//	id = m_neighborhoodSearch->add_static_point_set(x, (int)n);
			//}
			m_neighborhoodSearch->set_active_search(id, search_neighbors, find_neighbors);
			
			//// jfernandez: This is not necessary anymore. I leave it here while we debug the Akinci boundary handling
			//if (!is_dynamic)
			//	m_neighborhoodSearch->set_active_search((int) id, (int) id, false);
#else		
			id = m_neighborhoodSearch->add_point_set(x, n, is_dynamic, search_neighbors, find_neighbors);
#endif	

			if (id >= m_mapPointSet_UserData.size())
				m_mapPointSet_UserData.resize(id + 1);
			m_mapPointSet_UserData[id] = user_data;

			return id;
		}

#if defined(USE_TreeNSearch)
		unsigned int addPointSet(Real* x, const Real* radii, std::size_t n, bool is_dynamic, bool search_neighbors, bool find_neighbors, void* user_data)
		{
			unsigned int id = 0;
			if (is_dynamic) {
				id = m_neighborhoodSearch->add_point_set(x, radii, (int)n);
			}
			else 
			{
				//std::cout << "TreeNSearch warning: Cannot add static sets with variable search radius because static sets cannot search into other sets. The static was added as points to be found." << std::endl;
				//id = m_neighborhoodSearch->add_static_point_set(x, (int)n);
				id = m_neighborhoodSearch->add_point_set(x, (int)n);
			}

			m_neighborhoodSearch->set_active_search(id, search_neighbors, find_neighbors);

			if (id >= m_mapPointSet_UserData.size())
				m_mapPointSet_UserData.resize(id + 1);
			m_mapPointSet_UserData[id] = user_data;

			return id;
		}
#endif	

		void resizeSet(unsigned int pointSetIndex, Real const* x, std::size_t n)
		{
#if defined(USE_TreeNSearch)
			m_neighborhoodSearch->resize_point_set((int)pointSetIndex, x, (int)n);
#else
			m_neighborhoodSearch->resize_point_set(pointSetIndex, x, n);
#endif
		}

		void reset()
		{
#if defined(USE_CompactNSearch)
			m_neighborhoodSearch->reset();
#endif
		}

		void setSearchRadius(const Real radius)
		{
#if defined(USE_cuNSearch)
			m_neighborhoodSearch->set_radius(radius);
#elif defined(USE_CompactNSearch)
			m_neighborhoodSearch->set_radius(radius);
#elif defined(USE_TreeNSearch)
			m_neighborhoodSearch->set_search_radius(radius);
#endif
		}

		FORCE_INLINE unsigned int numberOfPointsInSet(const unsigned int pointSetIndex) const
		{
#if defined(USE_TreeNSearch)
			return static_cast<unsigned int>(m_neighborhoodSearch->get_n_points_in_set(pointSetIndex));
#else
			return static_cast<unsigned int>(m_neighborhoodSearch->point_set(pointSetIndex).n_points());
#endif
		}

		FORCE_INLINE unsigned int numberOfPointSets() const
		{
#if defined(USE_TreeNSearch)
			return static_cast<unsigned int>(m_neighborhoodSearch->get_n_sets());
#else
			return static_cast<unsigned int>(m_neighborhoodSearch->n_point_sets());
#endif
		}

		FORCE_INLINE unsigned int numberOfNeighbors(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index) const
		{
#if defined(USE_TreeNSearch)
			return static_cast<unsigned int>(m_neighborhoodSearch->get_neighborlist(pointSetIndex, neighborPointSetIndex, index).size());
#else
			return static_cast<unsigned int>(m_neighborhoodSearch->point_set(pointSetIndex).n_neighbors(neighborPointSetIndex, index));
#endif
		}

		FORCE_INLINE unsigned int getNeighbor(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index, const unsigned int k) const
		{
#if defined(USE_TreeNSearch)
			return m_neighborhoodSearch->get_neighborlist(pointSetIndex, neighborPointSetIndex, index)[k];
#else
			return m_neighborhoodSearch->point_set(pointSetIndex).neighbor(neighborPointSetIndex, index, k);
#endif
		}

		FORCE_INLINE const unsigned int* getNeighborList(const unsigned int pointSetIndex, const unsigned int neighborPointSetIndex, const unsigned int index) const
		{
#if defined(USE_cuNSearch)
			return m_neighborhoodSearch->point_set(pointSetIndex).neighbor_list(neighborPointSetIndex, index);
#elif defined(USE_CompactNSearch)
			return m_neighborhoodSearch->point_set(pointSetIndex).neighbor_list(neighborPointSetIndex, index).data();
#elif defined(USE_TreeNSearch)
			return reinterpret_cast<const unsigned int*>(m_neighborhoodSearch->get_neighborlist(pointSetIndex, neighborPointSetIndex, index).get_ptr());
			//return m_neighborhoodSearch->get_neighbor_list_ptr(pointSetIndex, neighborPointSetIndex, index);
#endif
		}
	};
}