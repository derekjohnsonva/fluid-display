<script>
	import { onMount } from 'svelte';
	import { invoke } from '@tauri-apps/api/tauri';

	// https://simeydotme.github.io/svelte-range-slider-pips/
	import RangeSlider from 'svelte-range-slider-pips';
	// Fluid Data
	let heat_val = [5];
	let heat_range = [0, 100];
	let bouyancy_val = [200];
	let bouyancy_range = [0, 500];
	let viscosity_range = [0, 500];
	let visc_val = [0];
	let diff_range = [0, 500];
	let diff_val = [0];
	// Screen Data
	let fluid_width = 100;
	let fluid_height = 100;

	let canvas;
	const invoke_func = window.__TAURI__.invoke;

	// draw a border around the canvas
	onMount(onload);
	async function onload() {
		// initialize the staggered grid
		await invoke_func('initialize_state', {
			width: fluid_width,
			height: fluid_height
		});
		const ctx = canvas.getContext('2d');
		let frame = setTimeout(loop);

		async function loop(t) {
			let data = await invoke_func('get_grid_colors', {})
			let data_array = Uint8ClampedArray.from(data);
			const imageData = new ImageData(data_array, fluid_width, fluid_height);
			ctx.putImageData(imageData, 0, 0);
			let result = await invoke_func('take_step', {})
			frame = setTimeout(loop);
		}
		return () => {
			clearTimeout(frame);
		};
	}
</script>

<div class="screen">
	<div class="form">
		<div class="inputs">
			<div class="slider">
				<p>
					Heat Transfer Rate: {heat_val / 100}
				</p>
				<RangeSlider
					bind:values={heat_val}
					min={heat_range[0]}
					max={heat_range[1]}
					on:stop={() => {
						invoke_func('update_heat_transfer', { newVal: heat_val[0] });
						console.debug('updated heat transfer rate');
					}}
				/>
			</div>
			<div class="slider">
				<p>
					Bouyancy: {bouyancy_val / 100}
				</p>
				<RangeSlider bind:values={bouyancy_val} min={bouyancy_range[0]} max={bouyancy_range[1]}
					on:stop={() => {
						invoke_func('update_bouyancy', { newVal: bouyancy_val[0] });
						console.debug('updated bouyancy');
					}}
				/>
			</div>
			<div class="slider">
				<p>
					Viscosity: {visc_val / 100}
				</p>
				<RangeSlider bind:values={visc_val} min={viscosity_range[0]} max={viscosity_range[1]} 
				on:stop={() => {
					invoke_func('update_viscosity', { newVal: visc_val[0] });
					console.debug('updated viscosity');
				}}
				/>
			</div>
			<div class="slider">
				<p>
					Diffusion: {diff_val / 100}
				</p>
				<RangeSlider bind:values={diff_val} min={diff_range[0]} max={diff_range[1]}
				on:stop={() => {
					invoke_func('update_diffusion', { newVal: diff_val[0] });
					console.debug('updated diffusion');
				}}
				/>
			</div>
		</div>
	</div>
	<div class="fluid">
		<div class="center">
			<canvas bind:this={canvas} height={fluid_height} width={fluid_width} />
		</div>
	</div>
</div>

<style>
	.screen {
		display: flex;
		flex-direction: row;
		justify-content: center;
		height: 100vh;
	}
	.form {
		margin: 0 auto;
		flex: 1;
		align-items: center;
		flex-direction: column;
		justify-content: center;
	}
	.inputs {
		justify-content: center;
		align-items: center;
		max-width: 16rem;
	}
	.fluid {
		flex: 1;
		background-color: black;
		justify-content: center;
		align-items: center;
		display: flex;
	}
	.slider {
		margin: 10px;
		height: 100%;
	}
</style>
